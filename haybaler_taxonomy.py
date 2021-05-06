# script to add taxonomy data to haybaler output .csv using pytaxonkit
# Sophia Poertner, Jan - May 2021
# Usage: conda activate haybaler
# bash run_haybaler_tax.sh
# Usage: python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /mnt/ngsnfs/seqres/metagenref/ --test_reference True
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


def find_species(csv, input_file, input_path):
    # try to filter the species out of the reference name
    species = []
    for organism in csv.index:
        organism = organism.replace("organism_", "")
        organism = organism.replace("1_1_1", "Homo_sapiens_")
        organism = re.sub(r'^.*?:', '', organism)  # replace everything before an ":" with nothing
        # find: two words containing just letters separated with "_". The first word starts with a capital letter
        m = re.findall('[A-Z][a-z]+_[a-z]+', organism)
        if len(m) != 0:
            m[0] = m[0].replace("_", " ")
            name = (m[0])
        else:
            name = "NOT KNOWN"
        species.append(name)
    return species


def find_genus(csv):
    genus = []
    for organism in csv.index:
        organism = organism.replace("organism_", "")
        organism = re.sub(r'^.*?:', '', organism)  # replace everything before an ":" with nothing
        organism = organism.replace("1_1_1", "Homo_sapiens_")
        split = organism.split(sep="_")
        if split[0] in ("NC", "AC", "NZ", "ENA"):
            del split[0]
        for element in split:
            if not has_numbers(element) and element:
                if not re.search("[Hh]uman", element):
                    name = element
                    break
                else:  # for e.g. NC_001352_1_Human_papillomavirus___2__complete_genome_VIR
                    next_element_index = split.index(element) + 1
                    next_element = split[next_element_index]
                    name = element + " " + next_element
                    break
            else:
                name = "NOT KNOWN"
        genus.append(name)
    return genus


def find_double_taxid(df, taxonomy_name):
    # In the package pytaxonkit is a phenomenon which makes it hard to assign a lineage to the original reference name
    # The method name2taxid assigns a TaxID to every input name. But a few different organisms share the same name so
    # a few names there have 2 or more TaxIDs. name2taxid outputs all possible names so sometimes there are 2 or 3
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
            # if three following organisms have the same Name but different TaxIDs
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
    new_df_2 = pd.DataFrame(columns=columns)
    # now filter for same input name but different output name
    # always take the organism with the same output and input name
    for index, row in new_df.iterrows():
        # skip this  iteration
        if double_taxid:
            double_taxid = False
            continue
        if index + 1 < len(new_df):
            # if two follwing organisms have the same input name but different output names
            if row["name2taxid_name"] == new_df["name2taxid_name"][index + 1] and row["Name"] != new_df["Name"][index + 1]:
                double_taxid = True
                # check if one of the output names has the same name as the input name, if yes add this organism to the
                # new df, else add just the name and leave the lineage empty
                if row["Name"] == row["name2taxid_name"]:
                    new_df_2 = new_df_2.append(row)
                elif new_df["Name"][index + 1] == row["name2taxid_name"]:
                    new_df_2 = new_df_2.append(new_df.iloc[index + 1, :])
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df_2 = pd.concat([new_df_2, one_row])
            else:
                new_df_2 = new_df_2.append(row)
        else:
            new_df_2 = new_df_2.append(row)
    new_df_2 = new_df_2.reset_index(drop=True)
    return new_df_2


def get_taxonomy(rank):
    taxonomy = pytaxonkit.name2taxid(rank)
    taxids = taxonomy["TaxID"].to_list()
    lineage = pytaxonkit.lineage(taxids)
    filtered_csv = find_double_taxid(lineage, taxonomy["Name"])
    return filtered_csv, taxonomy


def save_csv(csv, path, name):
    if "haybaler" in name:
        csv.to_csv(path + "/" + name.replace("haybaler", "haybaler_taxa"), sep="\t")
    else:
        sys.exit("ERROR: Input file {} has an incompatible file name. Needs a *haybaler.csv as input otherwise the "
                 "inputfile gets overwritten.".format(name))


def report(taxa, filtered_csv, taxonomy, csv, input_file, rank):
    if len(taxa) == len(filtered_csv):
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
        # print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit

        # insert taxa and lineage as series in the df and replace nan with ";;;;;;
        genus_series = pd.Series(taxa)
        lineage_series = pd.Series(filtered_csv["Lineage"].values)
        csv.insert(loc=0, column=rank, value=genus_series.values)
        csv.insert(loc=1, column='lineage_' + rank, value=lineage_series.values)
        csv["lineage_" + rank] = csv["lineage_" + rank].fillna(";;;;;;")

        total_chr = len(taxa)
        chr_not_work = len(taxonomy[taxonomy["TaxID"].isna()])
        chr_work = total_chr - chr_not_work
        chr_no_lineage = csv.loc[csv["lineage_" + rank] == ";;;;;;"].shape[0]
        chr_lineage = total_chr - chr_no_lineage
        print("reference tested:", input_file)
        print(total_chr, "total chromosomes,", rank, "found for", chr_work, "chromosomes, no", rank,  "found for",
              chr_not_work, "chromosomes")
        print("Found", rank, " for", chr_work / total_chr, "of the reference")
        print("The input and the output for lineage are the same. Filtering has probably worked.")
        print(total_chr, "total chromosomes,", "lineage found for", chr_lineage, "chromosomes, no lineage found for",
              chr_no_lineage, "chromosomes")
        print("Found lineage for", chr_lineage / total_chr, "of the reference")
        print("")
    else:
        print("The input for lineage is", len(taxa), "organisms long. The output is", len(filtered_csv),
              "long. In the process of getting taxid, lineage and filtering something must have gone wrong")
        print("")
        # uncomment the following lines to compare the lineages in different filter steps. Good for bugs
        # all_steps = pd.concat((pd.Series(csv.index.values), pd.Series(genus)), axis=1, ignore_index=True)
        # all_steps = pd.concat((all_steps, taxonomy[["Name", "TaxID"]]), axis=1, ignore_index=True)
        # all_steps = pd.concat((all_steps, lineage["Name"]), axis=1, ignore_index=True)
        # all_steps = pd.concat((all_steps, filtered_csv["Name"]), axis=1, ignore_index=True)
        # all_steps.columns = ["reference_name", "filtered_genus", "name_to_Taxid_Name", "name_to_taxid_TaxID", "lineage_Name", "filtered_csv_Name"]
        # all_steps.to_csv("compare_filter_steps.txt", index=False)
        # print(all_steps)


def add_taxonomy_to_df(csv, taxonomy, filtered_csv, df, path, file, rank):
    nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())
    # replace every name that produces as NAN output in pytaxonkit with "NOT KNOWN"
    genus_less_nan = ["NOT KNOWN" if name in nan_list else name for name in df]
    genus_series = pd.Series(genus_less_nan)
    lineage_series = pd.Series(filtered_csv["Lineage"].values)
    # insert genus and species name and lineage at the end of the csv
    if len(csv) == len(lineage_series):
        csv.insert(loc=len(csv.columns), column=rank + '_name', value=genus_series.values)
        csv.insert(loc=len(csv.columns), column=rank + '_lineage', value=lineage_series.values)
        save_csv(csv, path, file)
    else:
        print("The input for lineage is", len(df), "organisms long. The output is", len(filtered_csv),
              "long. In the process of getting taxid, lineage and filtering something must have gone wrong")
        print("")
    return csv


@click.command()
@click.option('--input_file', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Path of the input file, use . for current directory', required=True)
@click.option('--test_reference', '-t', help='Set to True if references should be tested. Default = False', default=False)
def main(input_file, input_path, test_reference):
    # pd.set_option('display.max_rows', 100000)
    if not test_reference:
        csv = pd.read_csv(input_path + "/" + input_file, sep="\t", index_col=0, header=0)
    else:
        csv = pd.read_csv(input_path + "/" + input_file, sep="\t", index_col=0, header=None)
    genus = find_genus(csv)
    filtered_genus, taxonomy_genus = get_taxonomy(genus)
    species = find_species(csv, input_file, input_path)
    filtered_species, taxonomy_species = get_taxonomy(species)
    if not test_reference:
        csv = add_taxonomy_to_df(csv, taxonomy_genus, filtered_genus, genus, input_path, input_file, rank="genus")
        add_taxonomy_to_df(csv, taxonomy_species, filtered_species, species, input_path, input_file, rank="species")
    else:
        report(genus, filtered_genus, taxonomy_genus, csv, input_file, rank="genus")
        report(species, filtered_species, taxonomy_species, csv, input_file, rank="species")


if __name__ == "__main__":
    main()
