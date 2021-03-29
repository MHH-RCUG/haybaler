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
def find_species(csv):
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
        print(name)
        species.append(name)
    taxonomy = pytaxonkit.name2taxid(species)
    nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
    genus_series = pd.Series(species)
    csv.insert(loc=0, column='genus', value=genus_series.values)
    print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit
    total_chr = len(species)
    chr_not_work = len(taxonomy[taxonomy["TaxID"].isna()])
    chr_work = total_chr - chr_not_work
    print("reference tested:", "test")
    print(total_chr, "total chromosomes,", chr_work, "chromosomes work,", chr_not_work, "do not work")
    print(chr_work / total_chr, "of the reference works,", chr_not_work / total_chr, "works not")
    print("")
    # print(species)
    # print(len(species))
    # print(species.count("NOT KNOWN"))
    return species


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


def save_csv(csv, path, name):
    if "haybaler" in name:
        csv.to_csv(path + "/" + name.replace("haybaler", "haybaler_genus"), sep="\t")
    else:
        sys.exit("ERROR: Input file {} has an incompatible file name. Needs a *haybaler.csv as input otherwise the inputfile "
                 "gets overwritten.".format(name))


@click.command()
@click.option('--input_file', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Path of the input file, use . for current directory', required=True)
def main(input_file, input_path):
    # Mode for testing References. True or False
    test_references = False
    pd.set_option('display.max_rows', 100000)
    csv = read_csv(input_file, input_path)
    # find_species(csv)  # work in progress
    genus = shorten_organism_names(csv)
    taxonomy = pytaxonkit.name2taxid(genus)
    if not test_references:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())
        # replace every name that produces as NAN output in pytaxonkit with "NOT KNOWN"
        genus_less_nan = ["NOT KNOWN" if name in nan_list else name for name in genus]
        genus_series = pd.Series(genus_less_nan)
        csv.insert(loc=0, column='genus', value=genus_series.values)
        save_csv(csv, input_path, input_file)
    else:
        nan_list = set(taxonomy[taxonomy["TaxID"].isna()]["Name"].to_list())  # list of names that produce NAN
        genus_series = pd.Series(genus)
        csv.insert(loc=0, column='genus', value=genus_series.values)
        print(csv[csv['genus'].isin(nan_list)]["genus"])  # print everything that didn't work with pytaxonkit
        # print(taxonomy[taxonomy["TaxID"].isna()])  # print everything that didn't worked with pytaxonkit (old)
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
