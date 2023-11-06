import pandas as pd
import subprocess
import re
import argparse
from collections import defaultdict

config = {
    (1, 26): "FR1",
    (27, 38): "CDR1",
    (39, 55): "FR2",
    (56, 65): "CDR2",
    (66, 104): "FR3",
    (105, 117): "CDR3",
    (118, 129): "FR4",
}


def extract_sequences(input_file):
    input_excel_file = pd.read_excel(input_file, sheet_name=1)
    input_vhsequences_data = input_excel_file.iloc[:, 3].tolist()
    input_vlsequences_data = input_excel_file.iloc[:, 4].tolist()

    with open('output_vh_file.fasta', 'w') as vh_file, open('output_vl_file.fasta', 'w') as vl_file:
        for i, sequence in enumerate(input_vhsequences_data):
            vh_file.write(f'>Sequence{i + 1}\n')
            vh_file.write(f'{sequence}\n')

        for i, sequence in enumerate(input_vlsequences_data):
            vl_file.write(f'>Sequence{i + 1}\n')
            vl_file.write(f'{sequence}\n')


def run_anarci(input_file, output_file):
    subprocess.run(["ANARCI", "-i", input_file, "-s", "imgt", "-o", output_file])


def add_annotation(input_file):
    aa_dict = defaultdict(lambda: defaultdict(str))
    num_of_seq = 0
    with open(input_file, 'r') as in_file:
        for line in in_file:
            if line.startswith('# Sequence'):
                num_of_seq += 1
            elif line.startswith(('L ', 'H ')):
                num = line.split()[1]
                length = len(line.split())
                for k, v in config.items():
                    if k[0] <= int(num) <= k[1]:
                        if length == 3:
                            aa_dict[f'Sequence{num_of_seq}'][v] += line.split()[2]
                        elif length == 4:
                            aa_dict[f'Sequence{num_of_seq}'][v] += line.split()[3]
                        else:
                            print('The content is not correct in the line, please check the input file!')

    seq = [f'Sequence{i}' for i in range(1, num_of_seq + 1)]
    data = defaultdict(list)
    data['Sequence'] = seq
    for k, v in aa_dict.items():
        for k1, v1 in v.items():
            data[k1].append(v1)
    return data


def annotate_sequences(input_file, output_file):
    data = add_annotation(input_file)
    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--input_file", help="输入文件")
    parser.add_argument("-out", "--output_file", help="输出文件")
    args = parser.parse_args()

    if args.input_file and args.output_file:
        extract_sequences(args.input_file)
        run_anarci("output_vh_file.fasta", "output_vh_numbered.txt")
        annotate_sequences("output_vh_numbered.txt", args.output_file)
    else:
        print("请输入输入文件和输出文件的路径")


if __name__ == "__main__":
    main()
