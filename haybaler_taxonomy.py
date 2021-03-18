# script to add taxaonomy data to habyaler output .csv useing pytaxonkit
# Sophia Poertner, Jan 2021

import pytaxonkit
import pandas as pd
import click
import re


def hasNumbers(inputString):
    return any(char.isdigit() for char in inputString)


def read_csv(file, path):
    csv = pd.read_csv(path + "/" + file, decimal=",", index_col=0, sep='\t')
    return csv


def shorten_organism_names(csv):
    genus = []
    for organism in csv.index:
        organism = organism.replace("1_1_1_", "chr")
        split = organism.split(sep="_")
        name = find_genus(split, organism)
        genus.append(name)
        print(organism, "    ", name)
    return genus


def find_genus(split, refseq_name):
    # for human chromosomes
    if re.search("^chr", split[0]):
        genus = split[0]
    else: 
        if split[0] in ("NC", "AC", "NZ", "ENA"):
            del split[0]
        for element in split:
            if not hasNumbers(element):
                if not re.search("[Hh]uman", element):
                    genus = element
                else: # for e.g. NC_001352_1_Human_papillomavirus___2__complete_genome_VIR
                    next_element_index = split.index(element) + 1
                    next_element = split[next_element_index]
                    genus = element + " " + next_element
                break
    if "genus" not in locals():
        print("It was not possible to detect the genus for ", refseq_name)
        genus = "UNKNOWN"
        print(genus)
    return genus


@click.command()
@click.option('--input_file', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Path of the input file', required=True)
def main(input_file, input_path):

    pd.set_option('display.max_rows', 100000)
    print(pytaxonkit.name2taxid(["Human papillomavirus", "Human herpesvirus"]))
    csv = read_csv(input_file, input_path)
    genus = shorten_organism_names(csv)
    taxonomy = pytaxonkit.name2taxid(genus)
    # taxonomy = pd.concat([taxonomy, pd.Series(genus)], axis=1, join="outer")
    print(taxonomy)
    print(taxonomy[taxonomy["TaxID"].isna()])
    result = pytaxonkit.lineage(taxonomy["TaxID"])
    print(result[['TaxID', 'Name', 'Lineage']])
    # taxonomy = pd.concat([taxonomy, pytaxonkit.lineage(taxonomy["TaxID"])], axis=1, join="outer")
    # print(taxonomy)
    # taxonomy.to_csv("taxonomy.csv")
    # pd.Series(genus).to_csv("genus.csv")

if __name__ == "__main__":
    main()
