# script to add taxonomy data to haybaler output .csv using pytaxonkit
# Sophia Poertner, Jan - March 2021
# Usage: conda activate haybaler
# Usage: python3 haybaler_taxonomy.py  -i 2021_02_human_bact_fungi_vir_masked.fa.fai -p /lager2/rcug/seqres/metagenref/
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
    result = pytaxonkit.lineage(taxids)
    # print(result)
    # print(len(csv), len(result["FullLineage"]), len(genus_series), len(taxids), len(taxonomy))
    # print(taxonomy, genus_series)
    # print(taxonomy)
    value = result["FullLineage"].values
    # value = result["Name"].values
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
        # print(organism, "    ", name)
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
                    # taxid = pytaxonkit.name2taxid([genus])
                    # print(taxid["Rank"][0])
                    # print(type(taxid["Rank"][0]))
                    # if taxid["Rank"][0].isna():
                    #     print("yes is na")
                    #     print(taxid["Rank"][0])
                    break
                else:  # for e.g. NC_001352_1_Human_papillomavirus___2__complete_genome_VIR
                    next_element_index = split.index(element) + 1
                    next_element = split[next_element_index]
                    genus = element + " " + next_element
                    # taxid = (pytaxonkit.name2taxid([genus]))
                    # if taxid["Rank"][0].isna():
                    #     print(taxid["Rank"][0])
                    # print(genus)
                    break
    if "genus" not in locals():
        # print("It was not possible to detect the genus for ", refseq_name)
        genus = "NOT KNOWN"
    return genus


def find_doubble_taxa(df, genus, taxonomy_name):
    # print(df)
    df["name2taxid_name"] = taxonomy_name
    # df.insert(5, "name2taxid_name", taxonomy_name)
    doubble_taxa = False
    doubble_taxa_2 = False
    columns = list(df.columns.values)
    new_df = pd.DataFrame(columns=columns)
    not_df = pd.DataFrame(columns=columns)
    for index, row in df.iterrows():
        if doubble_taxa_2:
            doubble_taxa_2 = False
            doubble_taxa = True
            continue
        if doubble_taxa:
            doubble_taxa = False
            continue
        if index + 1 != len(df):
            # elif all have the same Name but different TaxIDs
            if row["Name"] == df["Name"][index + 1] and row["Name"] == df["Name"][index + 2] and row["TaxID"] != df["TaxID"][index + 1] and row["TaxID"] != df["TaxID"][index + 2] and df["TaxID"][index + 2] != df["TaxID"][index + 1]:
                kingdome_1 = row["Lineage"].split(";")[0]
                kingdome_2 = df["Lineage"][index + 1].split(";")[0]
                kingdome_3 = df["Lineage"][index + 1].split(";")[0]
                doubble_taxa_2 = True
                if kingdome_1 == "Bacteria" and kingdome_2 != "Bacteria" and kingdome_3 != "Bacteria":
                    new_df = new_df.append(row)
                    not_df = not_df.append(df.iloc[index + 1, :])
                    not_df = not_df.append(df.iloc[index + 2, :])
                elif kingdome_2 == "Bacteria" and kingdome_1 != "Bacteria" and kingdome_3 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 1, :])
                    not_df = not_df.append(row)
                    not_df = not_df.append(df.iloc[index + 2, :])
                elif kingdome_3 == "Bacteria" and kingdome_1 != "Bacteria" and kingdome_2 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 2, :])
                    not_df = not_df.append(row)
                    not_df = not_df.append(df.iloc[index + 1, :])
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df = pd.concat([new_df, one_row])
                    not_df = not_df.append(row)
                    not_df = not_df.append(df.iloc[index + 1, :])
                    not_df = not_df.append(df.iloc[index + 2, :])
            elif row["Name"] == df["Name"][index + 1] and row["TaxID"] != df["TaxID"][index + 1]:
                # print(row["Name"], row["TaxID"], df["TaxID"][index + 1], row["Lineage"], df["Lineage"][index + 1])
                kingdome_1 = row["Lineage"].split(";")[0]
                kingdome_2 = df["Lineage"][index + 1].split(";")[0]
                doubble_taxa = True
                if kingdome_1 == "Bacteria" and kingdome_2 != "Bacteria":
                    new_df = new_df.append(row)
                    not_df = not_df.append(df.iloc[index + 1, :])
                elif kingdome_2 == "Bacteria" and kingdome_1 != "Bacteria":
                    new_df = new_df.append(df.iloc[index + 1, :])
                    not_df = not_df.append(row)
                else:
                    one_row = pd.DataFrame([row["Name"]], columns=["Name"])
                    new_df = pd.concat([new_df, one_row])
                    not_df = not_df.append(row)
                    not_df = not_df.append(df.iloc[index + 1, :])
            else:
                new_df = new_df.append(row)
        else:
            new_df = new_df.append(row)
    not_ids = set(not_df["TaxID"].values.tolist())
    new_df = new_df[~new_df['TaxID'].isin(not_ids)]
    new_df = new_df.reset_index(drop=True)
    new_df["Name"] = new_df["Name"].astype(str)
    columns = list(new_df.columns.values)
    not_df_new = pd.DataFrame(columns=columns)
    # print(new_df.loc[new_df['name2taxid_name'] == "Candidatus"])
    # print(new_df.iloc[1876]['Name'])
    # print(new_df.iloc[[1876]]["Name"])
    # print(type(new_df.iloc[1876]["Name"]))
    for index, row in new_df.iterrows():
        if index + 1 != len(new_df):
            if row["name2taxid_name"] == new_df["name2taxid_name"][index + 1] and row["Name"] != new_df["Name"][index + 1]:
                print(row["name2taxid_name"], row["Name"], new_df["Name"][index + 1])
                print(type(new_df["Name"][index + 1]), type(df["Name"][index + 1]), index)
                print("nan" == new_df["Name"][index + 1])
                if row["Name"] == row["name2taxid_name"]:
                    not_df_new = not_df_new.append(new_df.iloc[index + 1, :])
                if new_df["Name"][index + 1] == row["name2taxid_name"]:
                    not_df_new = not_df_new.append(row)
    not_ids_new = set(not_df_new["TaxID"].values.tolist())
    new_df = new_df[~new_df['TaxID'].isin(not_ids_new)]
    genus = pd.Series(genus)
    # print(not_ids)
    new_df_mini = pd.DataFrame({'genus': genus})
    new_df = new_df.reset_index(drop=True)
    # print(new_df["Lineage"])
    # new_df = pd.concat((new_df, genus), axis=1, ignore_index=True)
    # print(new_df[[3, 10, 2]])
    # print(new_df.columns.values)
    # print(new_df[["genus", "Name", "TaxID", "Lineage"]])
    # print(new_df)
    # print(len(new_df))
    # print(len(not_df))
    # print(len(df))
    # print("a", a)
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
def main(input_file, input_path):
    # Mode for testing References. True or False
    test_references = False
    pd.set_option('display.max_rows', 100000)
    csv = read_csv(input_file, input_path)
    # find_species(csv, input_file, input_path)  # work in progress
    genus = shorten_organism_names(csv)
    taxonomy = pytaxonkit.name2taxid(genus)
    taxids = taxonomy["TaxID"].to_list()
    result = pytaxonkit.lineage(taxids)
    # print(result)
    no_doubble_taxa = find_doubble_taxa(result, genus, taxonomy["Name"])  # [["TaxID", "Name", "Lineage"]]
    ###
    # print(pd.Series(csv.index.values))
    # print(type(pd.Series(genus)))
    # print(type(taxonomy[["Name", "TaxID"]]))
    # print(type(result["Name"]))
    # print(type(no_doubble_taxa["Name"]))
    all_steps = pd.concat((pd.Series(csv.index.values), pd.Series(genus)), axis=1, ignore_index=True)
    all_steps = pd.concat((all_steps, taxonomy[["Name", "TaxID"]]), axis=1, ignore_index=True)
    all_steps = pd.concat((all_steps, result["Name"]), axis=1, ignore_index=True)
    all_steps = pd.concat((all_steps, no_doubble_taxa["Name"]), axis=1, ignore_index=True)
    all_steps.columns = ["reference_name", "filtered_genus", "name_to_Taxid_Name", "name_to_taxid_TaxID", "lineage_Name", "no_doubble_taxa_Name"]
    all_steps.to_csv("compare_filter_steps.txt", index=False)
    # print(all_steps)
    # ###
    # print(len(genus))
    # print(result)
    if not test_references:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())
        # replace every name that produces as NAN output in pytaxonkit with "NOT KNOWN"
        genus_less_nan = ["NOT KNOWN" if name in nan_list else name for name in genus]
        genus_series = pd.Series(genus_less_nan)
        lineage_series = pd.Series(no_doubble_taxa["Lineage"].values)
        csv.insert(loc=0, column='genus', value=genus_series.values)
        csv.insert(loc=1, column='lineage', value=lineage_series.values)
        print(csv["lineage"])
        save_csv(csv, input_path, input_file)
    else:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
        genus_series = pd.Series(genus)
        # print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit
        # print(taxonomy[taxonomy["TaxID"].isna()])  # print everything that didn't worked with pytaxonkit (old)
        # print(taxonomy)
        total_chr = len(genus)
        chr_not_work = len(taxonomy[taxonomy["TaxID"].isna()])
        chr_work = total_chr - chr_not_work
        print("reference tested:", input_file)
        print(total_chr, "total chromosomes,", chr_work, "chromosomes were OK,", chr_not_work, "Did not work")
        print(chr_work / total_chr, "of the reference works,", chr_not_work / total_chr, "Did not work")
        print("")
    # print(result[['TaxID', 'Name', 'Lineage']])


if __name__ == "__main__":
    main()
