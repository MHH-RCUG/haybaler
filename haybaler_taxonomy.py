# script to add taxonomy data to haybaler output .csv using pytaxonkit
# Sophia Poertner, Jan - March 2021
# Usage: conda activate haybaler
# Usage: python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /mnt/ngsnfs/seqres/metagenref/
# Usage: python3 haybaler_taxonomy.py  -i RPMM_haybaler.csv -p ./control_dataset/haybaler_output/

import pytaxonkit
import pandas as pd
import click
import re
import sys


def has_numbers(input_string):
    return any(char.isdigit() for char in input_string)


def read_csv(file, path):
    csv = pd.read_csv(path + "/" + file, decimal=",", index_col=0, sep='\t')
    return csv


# work in progress - not used yet
def find_species(csv, input_file, input_path):
    species = []
    for organism in csv.index:
        organism = organism.replace("organism_", "")
        organism = organism.replace("1_1_1", "Homo_sapiens_")
        organism = re.sub(r'^.*?:', '', organism)
        m = re.findall('[A-Z][a-z]+_[a-z]+', organism)
        if len(m) != 0:
            m[0] = m[0].replace("_", " ")
            name = (m[0])
        else:
            name = "NOT KNOWN"
        # print(name)
        species.append(name)
    taxonomy = pytaxonkit.name2taxid(species)
    nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
    genus_series = pd.Series(species)
    csv.insert(loc=0, column='genus', value=genus_series.values)
    # print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit
    total_chr = len(species)
    chr_not_work = len(taxonomy[taxonomy["TaxID"].isna()])
    chr_work = total_chr - chr_not_work
    print("reference tested:", "test")
    print(total_chr, "total chromosomes,", chr_work, "chromosomes work,", chr_not_work, "do not work")
    print(chr_work / total_chr, "of the reference works,", chr_not_work / total_chr, "works not")
    print("")
    taxids = taxonomy["TaxID"].to_list()
    lineage = pytaxonkit.lineage(taxids)
    # print(lineage)
    # print(len(csv), len(lineage["FullLineage"]), len(genus_series), len(taxids), len(taxonomy))
    # print(taxonomy, genus_series)
    # print(taxonomy)
    value = lineage["FullLineage"].values
    # value = lineage["Name"].values
    csv.insert(loc=1, column="lineage", value=value)
    # print(csv)
    save_csv(csv, input_path, input_file.replace(".csv", "_species.csv"))
    # print(species)
    # print(species.count("NOT KNOWN"))


def shorten_organism_names(csv):
    genus = []
    for organism in csv.index:
        organism = organism.replace("organism_", "")
        organism = re.sub(r'^.*?:', '', organism)  # replace everything before an ":" with nothing
        split = organism.split(sep="_")
        name = find_genus(split, organism)
        genus.append(name)
    return genus


def find_genus(split, refseq_name):
    # for human chromosomes
    if re.search("^1_1_1_", refseq_name):
        genus = "Homo"
    else:
        if split[0] in ("NC", "AC", "NZ", "ENA"):
            del split[0]
        for element in split:
            if not has_numbers(element) and element:
                if not re.search("[Hh]uman", element):
                    genus = element
                    break
                else:  # for e.g. NC_001352_1_Human_papillomavirus___2__complete_genome_VIR
                    next_element_index = split.index(element) + 1
                    next_element = split[next_element_index]
                    genus = element + " " + next_element
                    break
    if "genus" not in locals():
        # print("It was not possible to detect the genus for ", refseq_name)
        genus = "NOT KNOWN"
    return genus


def find_double_taxid(df, taxonomy_name):
    # In the package pytaxonkit i a phenomenon which make it hard to assign a lineage to the original reference name
    # The method name2taxid assigns a TaxID to every input name. But a few different organisms share the same name so
    # a few names there are 2 or more TaxIDs. name2taxid outputs all possible names so sometimes there are 2 or 3
    # outputs for one input. We need to find out which TaxID is correct if possible.
    # -> filter for multiple lines with the same name but different TaxIDs
    # Additionally sometimes the method lineage outputs a different name than what was originally the input name so the
    # above filter does not work anymore. Eg: input: Bacteroidetes ouputs: Bacteroidetes; Bacteroidia
    # We also need the input names (taxonomy_name)
    # -> filter for multiple lines in the input with the same name but different names in the output
    df["name2taxid_name"] = taxonomy_name  # input name
    double_taxid = False  # True if one name has two TaxIDs
    triple_taxid = False  # True if one name has three TaxIDs
    columns = list(df.columns.values)
    new_df = pd.DataFrame(columns=columns)
    # first filter for same name but different taxID
    for index, row in df.iterrows():
        # skip this and the next iteration
        if triple_taxid:
            triple_taxid = False
            double_taxid = True
            continue
            # skip this  iteration
        if double_taxid:
            double_taxid = False
            continue
        if index + 2 < len(df):
            # if three have the same Name but different TaxIDs
            if row["Name"] == df["Name"][index + 1] and row["Name"] == df["Name"][index + 2] and row["TaxID"] != \
                    df["TaxID"][index + 1] and row["TaxID"] != df["TaxID"][index + 2] and df["TaxID"][index + 2] != \
                    df["TaxID"][index + 1]:
                domain_1 = row["Lineage"].split(";")[0]
                domain_2 = df["Lineage"][index + 1].split(";")[0]
                domain_3 = df["Lineage"][index + 1].split(";")[0]
                triple_taxid = True  # skip the next 2 iterations
                # check if only one of the organisms is a bacteria. If yes, add this organism to the new df, else add
                # just the name and leave the lineage empty
                if domain_1 == "Bacteria" and domain_2 != "Bacteria" and domain_3 != "Bacteria":
                    new_df = new_df.append(row)
                elif domain_2 == "Bacteria" and domain_1 != "Bacteria" and domain_3 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 1, :])
                elif domain_3 == "Bacteria" and domain_1 != "Bacteria" and domain_2 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 2, :])
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df = pd.concat([new_df, one_row])
            # if two have the same Name but different TaxIDs
            elif row["Name"] == df["Name"][index + 1] and row["TaxID"] != df["TaxID"][index + 1]:
                domain_1 = row["Lineage"].split(";")[0]
                domain_2 = df["Lineage"][index + 1].split(";")[0]
                double_taxid = True  # skip next iteration
                # check if only one of the organisms is a bacteria. If yes, add this organism to the new df, else add
                # just the name and leave the lineage empty
                if domain_1 == "Bacteria" and domain_2 != "Bacteria":
                    new_df = new_df.append(row)
                elif domain_2 == "Bacteria" and domain_1 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 1, :])
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df = pd.concat([new_df, one_row])
            else:
                new_df = new_df.append(row)
        elif index + 1 < len(df):
            # if two have the same Name but different TaxIDs
            if row["Name"] == df["Name"][index + 1] and row["TaxID"] != df["TaxID"][index + 1]:
                domain_1 = row["Lineage"].split(";")[0]
                domain_2 = df["Lineage"][index + 1].split(";")[0]
                double_taxid = True  # skip next iteration
                # check if only one of the organisms is a bacteria. If yes, add this organism to the new df, else add
                # just the name and leave the lineage empty
                if domain_1 == "Bacteria" and domain_2 != "Bacteria":
                    new_df = new_df.append(row)
                elif domain_2 == "Bacteria" and domain_1 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 1, :])
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df = pd.concat([new_df, one_row])
            else:
                new_df = new_df.append(row)
        else:
            new_df = new_df.append(row)
    new_df = new_df.reset_index(drop=True)
    new_df["Name"] = new_df["Name"].astype(str)
    columns = list(new_df.columns.values)
    df_wrong_organisms = pd.DataFrame(columns=columns)  # df for wrong organisms
    # now filter for same input name but different output name
    # always take the organism with the same output and input name
    for index, row in new_df.iterrows():
        if index + 1 != len(new_df):
            if row["name2taxid_name"] == new_df["name2taxid_name"][index + 1] and row["Name"] != new_df["Name"][index + 1]:
                if row["Name"] == row["name2taxid_name"]:
                    df_wrong_organisms = df_wrong_organisms.append(new_df.iloc[index + 1, :])  # add wrong organisms to df_wrong_organisms
                if new_df["Name"][index + 1] == row["name2taxid_name"]:
                    df_wrong_organisms = df_wrong_organisms.append(row)  # add wrong organisms to the df_wrong_organisms
    wrong_ids = set(df_wrong_organisms["TaxID"].values.tolist())  # get the IDs for the wrong organisms
    new_df = new_df[~new_df['TaxID'].isin(wrong_ids)]  # new df which does not contain the wrong organisms
    new_df = new_df.reset_index(drop=True)
    return new_df


def save_csv(csv, path, name):
    if "haybaler" in name:
        csv.to_csv(path + "/" + name.replace("haybaler", "haybaler_genus"), sep="\t")
    else:
        sys.exit("ERROR: Input file {} has an incompatible file name. Needs a *haybaler.csv as input otherwise the "
                 "inputfile gets overwritten.".format(name))


@click.command()
@click.option('--input_file', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Path of the input file, use . for current directory', required=True)
@click.option('--test_reference', '-t', help='Set to True if references should be tested. Default = False', default=False)
def main(input_file, input_path, test_reference):
    # Mode for testing References. True or False
    # test_references = True
    pd.set_option('display.max_rows', 100000)
    csv = pd.read_csv(input_path + "/" + input_file, sep="\t", index_col=0, header=None)
    # find_species(csv, input_file, input_path)  # work in progress
    genus = shorten_organism_names(csv)
    taxonomy = pytaxonkit.name2taxid(genus)
    taxids = taxonomy["TaxID"].to_list()
    lineage = pytaxonkit.lineage(taxids)
    filtered_csv = find_double_taxid(lineage, taxonomy["Name"])  # [["TaxID", "Name", "Lineage"]]
    if not test_reference:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())
        # replace every name that produces as NAN output in pytaxonkit with "NOT KNOWN"
        genus_less_nan = ["NOT KNOWN" if name in nan_list else name for name in genus]
        genus_series = pd.Series(genus_less_nan)
        lineage_series = pd.Series(filtered_csv["Lineage"].values)
        csv.insert(loc=0, column='genus', value=genus_series.values)
        csv.insert(loc=1, column='lineage', value=lineage_series.values)
        save_csv(csv, input_path, input_file)
    else:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
        genus_series = pd.Series(genus)
        lineage_series = pd.Series(filtered_csv["Lineage"].values)
        csv.insert(loc=0, column='genus', value=genus_series.values)
        csv.insert(loc=1, column='lineage', value=lineage_series.values)
        csv["lineage"] = csv["lineage"].fillna(";;;;;;")
        # print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit
        total_chr = len(genus)
        chr_not_work = len(taxonomy[taxonomy["TaxID"].isna()])
        chr_work = total_chr - chr_not_work
        chr_no_lineage = csv.loc[csv["lineage"] == ";;;;;;"].shape[0]
        chr_lineage = total_chr - chr_no_lineage
        print("reference tested:", input_file)
        print(total_chr, "total chromosomes,", "genus found for", chr_work, "chromosomes, no genus found for",
              chr_not_work, "chromosomes")
        print("Found genus for", chr_work / total_chr, "of the reference")
        if len(genus) == len(filtered_csv):
            print("The input and the output for lineage are the same. Filtering has probably worked.")
            print(total_chr, "total chromosomes,", "lineage found for", chr_lineage, "chromosomes, no lineage found for",
                  chr_no_lineage, "chromosomes")
            print("Found lineage for", chr_lineage / total_chr, "of the reference")
            print("")
        else:
            print("The input for lineage is", len(genus), "organisms long. The output is", len(filtered_csv),
                  "long. In the process of getting taxid, lineage and filtering something must have gone wrong")
            print("")
            # compare the lineages in the different filter steps. Good for bug and reference testing
            # all_steps = pd.concat((pd.Series(csv.index.values), pd.Series(genus)), axis=1, ignore_index=True)
            # all_steps = pd.concat((all_steps, taxonomy[["Name", "TaxID"]]), axis=1, ignore_index=True)
            # all_steps = pd.concat((all_steps, lineage["Name"]), axis=1, ignore_index=True)
            # all_steps = pd.concat((all_steps, filtered_csv["Name"]), axis=1, ignore_index=True)
            # all_steps.columns = ["reference_name", "filtered_genus", "name_to_Taxid_Name", "name_to_taxid_TaxID", "lineage_Name", "filtered_csv_Name"]
            # all_steps.to_csv("compare_filter_steps.txt", index=False)
            # print(all_steps)


if __name__ == "__main__":
    main()
