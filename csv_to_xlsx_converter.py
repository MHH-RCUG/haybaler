import pandas
import numpy
import sys

# Fabian Charly Friedrich, April 2020
# Use the runbatch_csv_to_xlsx.sh script
# Or python3 csv_to_xlsx_converter.py <csv file>


def read_annot(file):
    lines = open(file, 'r').readlines()
    data = []

    for line in lines:
        split_line = line.strip().split(' ')

        tmp = []
        for i in range(3):
            try:
                tmp.append(split_line[i])
            except IndexError:
                tmp.append("")

        try:
            tmp.append(' '.join(split_line[3:]))
        except IndexError:
            tmp.append("")

        data.append(tmp)

    return pandas.DataFrame(numpy.array(data))


def main():
    file = sys.argv[1]
    df = read_annot(file)
    df.to_excel(file.replace('.csv', '.xlsx'), index=None, header=False)


if __name__ == '__main__':
    main()
