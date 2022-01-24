# Haybaler
# Sophia Poertner, Nov 2020 - April 2021

# Combine your Wochenende .bam.txt or reporting output from multiple samples into one matrix per stat.
# Usage: bash run_haybaler.sh


import pandas as pd
import click
import os
import re

version = "0.30 - April 2021"


# changelog
# 0.30 read all samples in one call. Filter out taxa with values below a readcount and RPMM limit
# 0.23 improve file input and arg handling
# 0.22 bugfix, correct gc_ref and chr_length for new chromosomes
# 0.21 fix ordering problems
# 0.20 add find_order and sort_new functions, so taxa with highest readcounts come first
# 0.11 add heatmap prep and R scripts
# 0.10 initial commits, improvements, testing


def read_csv(filename, filepath):
    file = pd.read_csv(filepath + '/' + filename, decimal=",", index_col=0)
    return file


def txt_to_df(filename, filepath):
    with open(filepath + '/' + filename) as infile, open('tmp.csv', 'w') as outfile:
        # add column names (not given in txt.file), save new file as temp outfile
        outfile.write("species,chr_length,read_count,unmapped_read_segments\n")
        # replace tabs with comma(tab separated to comma separated)
        for line in infile:
            outfile.write(" ".join(line.split()).replace(' ', ','))
            outfile.write("\n")
    file = pd.read_csv("tmp.csv", decimal=",", index_col=0)
    if os.path.exists("tmp.csv"):  # del tmp file outfile
        os.remove("tmp.csv")
    del file['unmapped_read_segments']  # unneeded column?
    return file


def join_dfs(file, name, path, column, input_name):
    sample = (input_name[:input_name.find(".")])  # shorten sample name
    sub_df = file[[column]].copy()  # new df with just the wanted column
    sub_df = sub_df.rename(columns={column: sample})  # rename column to sample name
    if os.path.isfile(path + "/" + column + "_" + name):  # if the file for the wanted stat already exists
        old = pd.read_csv(path + "/" + column + "_" + name, decimal=",", index_col=0, sep='\t')
        old.fillna(0.0, inplace=True)
        if sample not in old.columns:  # no double samples
            new_chr = []  # get chromosomes which are new in this sample
            for chromosome in file.index:
                if chromosome not in old.index:
                    new_chr.append(chromosome)
            # get a df with the chr_length and gc_ref from the new chromosomes
            if 'gc_ref' in file:
                new_chr_df = file.loc[new_chr, ['chr_length', 'gc_ref']]
            else:
                new_chr_df = file.loc[new_chr, ['chr_length']]
            old = old.append(new_chr_df)  # append the df with chr_length and gc_ref to the old df
            new = pd.concat([old, sub_df], axis=1, sort=False)  # add the new column to the old df
            if 'gc_ref' not in new and 'gc_ref' in file:
                gc = file[['gc_ref']].copy()
                new = pd.concat([new, gc], axis=1, sort=False)
                tmp = new['gc_ref'].to_list()
                del new['gc_ref']
                new.insert(1, 'gc_ref', tmp)
        else:
            new = old
    # if the file for the wanted stat does not exist, make this file containing the columns which are always the same
    # and the current sample
    else:
        if 'gc_ref' in file:
            new = file[['chr_length', 'gc_ref', column]].copy()
        else:
            new = file[['chr_length', column]].copy()
        new = new.rename(columns={column: sample})
    new.fillna(0.0, inplace=True)
    new = new.astype(float)
    new = new.round(2)
    return new


# calculate in which order the organisms should be in the output files.
# the organism with the most read count in all samples should come first
def find_order(df):
    samples = []  # list with all samples
    for column in df.columns:
        if column != 'chr_length' and column != 'gc_ref':
            samples.append(column)
    sum_organisms = []  # list of the sum form all samples for each organism (row sums)
    for organism in df.index:
        tmp_organism = []  # list of the stats from all samples for one organism
        for column in samples:
            tmp_organism.append(float(df.at[organism, column]))
        sum_organisms.append(sum(tmp_organism))
    df['sum_organisms'] = sum_organisms  # add a column with the sums to the df
    df = df.sort_values(by='sum_organisms', ascending=False)  # sort the df by the sums
    df = df.drop(['sum_organisms'], axis=1)  # delete the column with the sums
    order = df.index
    return df, order


# sort the new df so it fits the previous calculated order
def sort_new(df, order):
    order_df = pd.DataFrame(index=order)  # create an empty order_df with just the right orderer organisms as index
    new = pd.concat([order_df, df], axis=1, sort=False)  # concat the df on the order_df so it is in the right order too
    return new


def adding_species(path, column, name):
    # when concating two df's, the value species gets lost, so it needs to be added afterwards
    with open(path + "/" + column + "_" + name, 'r+') as f:
        content = f.read()
        if not content.split()[0] == "species":
            f.seek(0, 0)
            f.write(f"species" + content)


def get_taxa_to_exclude(file, limit, taxa_to_exclude, path):
    reason = file.replace("_haybaler.csv", "_below_" + str(limit))
    reason = reason.replace(path + "/", "")
    df = pd.read_csv(file, decimal=",", index_col=0, sep='\t')
    df = df.drop(['chr_length', 'gc_ref'], axis=1)
    # check which rows have all values below the limit, add them to the "taxa_to_exclude_list"
    taxa_to_exclude.extend(df[(df.astype('float') < limit).all(axis=1)].index)
    # create a df with the excluded taxa and why they are excluded
    taxa_to_exclude_index = df[(df.astype('float') < limit).all(axis=1)].index
    taxa_to_exclude_df = pd.DataFrame(index=taxa_to_exclude_index)
    taxa_to_exclude_df[reason] = "yes"
    return taxa_to_exclude, taxa_to_exclude_df


def exclude_taxa(file, path, taxa_to_exclude):
    df = pd.read_csv(path + "/" + file, decimal=",", index_col=0, sep='\t')  # read csv
    df = df[~df.index.isin(taxa_to_exclude)]  # exclude taxa
    df.to_csv(path + "/" + file, sep="\t")  # save again as csv


def shorten_names(output_path, col, output_file):
    short = pd.read_csv(output_path + "/" + col + "_" + output_file, index_col = 0, sep = "\t")
    rownames = list(short.index)
    for row in rownames:
        new_name = row
        split_name = row.split("_")
        for n in range(len(split_name)):
            if split_name[n] == "organism" or split_name[n] == "candidatus":
                new_name = split_name[n] + "_" + split_name[n+1] + "_" + split_name[n+2]
                new_name = subspecies(new_name, n, 3, split_name)
                short = change_name(new_name, row, short)
                break
            try:
                first_letter = ord(split_name[n][0])
                second_letter = ord(split_name[n][1])
                second_first = ord(split_name[n+1][0])
                if 64 < first_letter < 91 and 96 < second_letter < 123 and 96 < second_first < 123:
                    new_name = split_name[n] + "_" + split_name[n+1]
                    new_name = subspecies(new_name, n, 2, split_name)
                    short = change_name(new_name, row, short)
                    break
            except:
                pass
    save_name = output_path + "/" + col + "_" + output_file
    try:
        save_name = save_name.split(".")[-2] + "_short.csv"
    except:
        save_name = save_name + "_short.csv"
    index = short.index
    if index.is_unique:  # only saved if all row names are unique
        short.to_csv(save_name, sep='\t')


def subspecies(new_name, n, count,split_name):
    if split_name[n+count] == "subsp":
        add = split_name[n+count+1]
        while len(add) == 0:
            count += 1
            add = split_name[n+count+1]
        new_name = new_name + "_subsp_" + add
    return(new_name)


def change_name(new_name, row, short):
    short.rename(index={row:new_name}, inplace=True)
    return(short)


@click.command()
@click.option('--input_files', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Path of the input file', required=True)
@click.option('--output_path', '-op', help='Name of the output path')
@click.option('--output_file', '-o', help='Name of the output file')
@click.option('--readcount_limit', '-l', help='minimum amount of readcounts per sample. "Chromosomes with less than x '
                                              'reads in every sample are filtered out"! Default = 10', default=10)
@click.option('--rpmm_limit', '-r', help='minimum amount of RPMM per sample. "Chromosomes with less than x '
                                         'RPMM in every sample are filtered out"! Default = 300', default=300)
def main(input_files, input_path, output_path, output_file, readcount_limit, rpmm_limit):
    list_input_files = input_files.split(";")[1:]
    col_list = []
    # Debug prints messages on input and progress
    debug = False  # True or False
    if debug:
        print("INFO: Haybaler debug is on.")

    for input_file in list_input_files:
        try:
            if input_file.endswith('.csv'):
                if not re.search("rep.u*s.csv", input_file):
                    print("WARNING: Input file {0} does not match the typical file names. Only bam.txt, rep.s.csv and "
                          "rep.us.csv work as input files.".format(input_file))
                file = read_csv(input_file, input_path)
                if debug:
                    print(input_file)
            elif input_file.endswith('.txt'):
                if not re.search("bam.txt", input_file):
                    print("WARNING: Input file {0} does not match the typical file names. Only bam.txt, rep.s.csv and "
                          "rep.us.csv work as input files.".format(input_file))
                file = txt_to_df(input_file, input_path)
                if debug:
                    print(input_file)
            else:
                raise Exception(
                    "Inputfile {0} has the wrong file format. Only works for txt and csv".format(input_file))
        except FileNotFoundError:
            raise Exception("Failed to find or read input file: {0}".format(input_file))
        except AttributeError:
            raise Exception("No input file given. Please specify input file. Try --help for help")
        # make an own file for each stat, so for each column in the input file (or add the stats to existing files)
        for col in file.columns:
            if col != "chr_length" and col != "gc_ref":  # columns which are the same in every sample. Don't need extra file
                if debug:
                    print(col)
                df = join_dfs(file, output_file, output_path, col, input_file)
                if col == 'read_count':
                    df, order = find_order(df)
                elif 'order' in locals():
                    df = sort_new(df, order)
                else:
                    print(
                        "Sorting process for file {0} passed. It was not possible to calculate the correct order".format
                        (input_file))
                df.to_csv(output_path + "/" + col + "_" + output_file, sep='\t')
                adding_species(output_path, col, output_file)
                col_list.append(col)


    taxa_to_exclude = []
    excluded_taxa_readcount = None
    excluded_taxa_RPMM = None
    if os.path.isfile(output_path + "/" + "read_count_" + output_file):
        taxa_to_exclude, excluded_taxa_readcount = get_taxa_to_exclude(output_path + "/" + "read_count_" + output_file,
                                                                       readcount_limit, taxa_to_exclude, output_path)
    if os.path.isfile(output_path + "/" + "RPMM_" + output_file):
        taxa_to_exclude, excluded_taxa_RPMM = get_taxa_to_exclude(output_path + "/" + "RPMM_" + output_file, rpmm_limit,
                                                                  taxa_to_exclude, output_path)

    # check if the df of the excluded taxa from readcount and RPMM exist, concat them if both do
    if excluded_taxa_readcount is not None and excluded_taxa_RPMM is not None:
        excluded_taxa_df = pd.concat([excluded_taxa_readcount, excluded_taxa_RPMM], axis=1)
    elif excluded_taxa_readcount is not None and excluded_taxa_RPMM is None:
        excluded_taxa_df = excluded_taxa_readcount
    elif excluded_taxa_readcount is None and excluded_taxa_RPMM is not None:
        excluded_taxa_df = excluded_taxa_RPMM
    else:
        excluded_taxa_df = pd.DataFrame()
    excluded_taxa_df.fillna("no", inplace=True)
    excluded_taxa_df.to_csv(output_path + "/excluded_taxa.csv", sep="\t")  # save the df with the excluded taxa and the reason
    # when concating two df's, the value species gets lost, so it needs to be added afterwards
    with open(output_path + "/excluded_taxa.csv", 'r+') as f:
        content = f.read()
        if not content.split()[0] == "species":
            f.seek(0, 0)
            f.write(f"species" + content)
    for haybaler_csv in os.listdir(output_path):
        if haybaler_csv.endswith(output_file):
            exclude_taxa(haybaler_csv, output_path, taxa_to_exclude)  # exclude the taxa from the haybaler.csv

    # recreating all the output csv's with only species names as row names
    for col in list(set(col_list)):
        shorten_names(output_path, col, output_file)


if __name__ == '__main__':
    main()
