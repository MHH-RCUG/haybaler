# Haybaler
# Sophia Poertner, Nov 2020, Jan 2021

# Combine your Wochenende .bam.txt or reporting output from multiple samples into one matrix per stat.
# Usage: bash run_haybaler.sh


import pandas as pd
import click
import os.path
import re

version = "0.21 - Feb 2021"

# changelog
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
    sample = (input_name[:input_name.find(".")])   # shorten sample name
    sub_df = file[[column]].copy()  # new df with just the wanted column
    sub_df = sub_df.rename(columns={column: sample})  # rename column to sample name
    if os.path.isfile(path + "/" + column + "_" + name):  # if the file for the wanted stat already exists
        old = pd.read_csv(path + "/" + column + "_" + name, decimal=",", index_col=0, sep='\t')
        old.fillna(0.0, inplace=True)
        if sample not in old.columns:  # no double samples
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


@click.command()
@click.option('--input_file', '-i', help='Name of the input file')
@click.option('--input_path', '-p', help='Path of the input file')
@click.option('--output_path', '-op', help='Name of the output path')
@click.option('--output_file', '-o', help='Name of the output file')
def main(input_file, input_path, output_path, output_file):
    # Debug prints messages on input and progress
    debug = False  # True or False
    if debug:
        print("INFO: Haybaler debug is on.")
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
            raise Exception("Inputfile {0} has the wrong file format. Only works for txt and csv".format(input_file))
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
                print("Sorting process for file {0} passed. It was not possible to calculate the correct order".format
                      (input_file))
            df.to_csv(output_path + "/" + col + "_" + output_file, sep='\t')
            adding_species(output_path, col, output_file)


if __name__ == '__main__':
    main()
