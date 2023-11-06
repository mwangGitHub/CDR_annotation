import pandas as pd
import subprocess
import re
import argparse
from collections import defaultdict
'''
相比script2.py添加了PDB ID信息

优化函数，三步合并成一步，完成从原始excel中提取vh和vl序列信息，转化成对应的fasta文件，
然后对fasta文件调用ANARCI软件，将aa序列编号，
最后将编号后的vh和vl进行标注区分FR和CDR区，执行命令如下：
python script3.py -in selected_data.xlsx -out out3.xlsx
'''

config = {
    (1, 26): "FR1",
    (27, 38): "CDR1",
    (39, 55): "FR2",
    (56, 65): "CDR2",
    (66, 104): "FR3",
    (105, 117): "CDR3",
    (118, 129): "FR4",
}
#第一步从原始xlsx文件中提取vh和vl序列信息，保存成fasta格式
def extract_sequences(input_file):
    input_excel_file = pd.read_excel(input_file, sheet_name=1) # 提取第2个sheet中的抗体数据
    input_vhsequences_data = input_excel_file.iloc[:, 3].tolist() # 重链在第3列，轻链在第4列
    input_vlsequences_data = input_excel_file.iloc[:, 4].tolist()

    with open('output_vh.fasta', 'w') as vh_file, open('output_vl.fasta', 'w') as vl_file:
        for i, sequence in enumerate(input_vhsequences_data):
            vh_file.write(f'>Sequence{i + 1}\n')
            vh_file.write(f'{sequence}\n')

        for i, sequence in enumerate(input_vlsequences_data):
            vl_file.write(f'>Sequence{i + 1}\n')
            vl_file.write(f'{sequence}\n')
    
#第二步，调用ANARCI包，对fasta文件中的aa进行编号
def run_anarci(input_file, output_file):
    subprocess.run(["ANARCI", "-i", input_file, "-s", "imgt", "-o", output_file])
#添加注释函数
def add_annotation(input_file):
    #检查输入是Heavh chain or Light chain
    aa_dict_vh = defaultdict(lambda: defaultdict(str))#字典嵌套字典，即将内层字典赋值给外层字典中的value
    num_of_seq_vh = 0
    aa_dict_vl = defaultdict(lambda: defaultdict(str))
    num_of_seq_vl = 0
 
    with open(input_file, 'r') as in_file:
        for line in in_file:
            if line.startswith('# Sequence'):
                if 'vh' in str(input_file):
                    num_of_seq_vh += 1
                elif 'vl' in str(input_file):
                    num_of_seq_vl += 1
            elif line.startswith(('L ', 'H ')):
                num = line.split()[1]#提取编号信息，位于每行的第2位;line的格式为L 1       D
                length = len(line.split())#获取行元素个数
                for k, v in config.items():
                    if k[0] <= int(num) <= k[1]:
                        if length == 3:#判断是否存在rearrangement
                            if 'vl' in str(input_file): #区分vl or vh
                                aa_dict_vl[f'Sequence{num_of_seq_vl}'][v] += line.split()[2]
                            else:
                                aa_dict_vh[f'Sequence{num_of_seq_vh}'][v] += line.split()[2]
                        elif length == 4:
                            if 'vl' in str(input_file):
                                aa_dict_vl[f'Sequence{num_of_seq_vl}'][v] += line.split()[3]
                            else:
                                aa_dict_vh[f'Sequence{num_of_seq_vh}'][v] += line.split()[3]
                        else:
                            print('The content is not correct in the line, please check the input file!')

    if 'vh' in str(input_file):
        seq_vh = [f'Sequence{i}' for i in range(1, num_of_seq_vh + 1)]
        data = defaultdict(list)
        data['Sequence'] = seq_vh
        for k, v in aa_dict_vh.items():
            for k1, v1 in v.items():
                data[k1].append(v1)

        df_vh = pd.DataFrame(data)
        return df_vh
    else:
        seq_vl = [f'Sequence{i}' for i in range(1, num_of_seq_vl + 1)]
        data = defaultdict(list)
        data['Sequence'] = seq_vl
        for k, v in aa_dict_vl.items():
            for k1, v1 in v.items():
                data[k1].append(v1)
        df_vl = pd.DataFrame(data)
        return df_vl

def annotate_seq(input_file, output_file):
    #添加原始xlsx文件中的PDB ID信息到输出文件中
    input_excel_file = pd.read_excel(input_file, sheet_name=1)
    input_pdbid_data = input_excel_file.iloc[:, 1].tolist()

    run_anarci("output_vh.fasta", "out_numbered_vh.txt")
    run_anarci("output_vl.fasta", "out_numbered_vl.txt")    
    df_vh = add_annotation("out_numbered_vh.txt")
    df_vl = add_annotation("out_numbered_vl.txt")
    df_vh.insert(1, 'PDB ID', input_pdbid_data)
    df_vl.insert(1, 'PDB ID', input_pdbid_data)

    with pd.ExcelWriter(output_file) as writer:
        df_vh.to_excel(writer, sheet_name="Heavy_chain", index=False)
        df_vl.to_excel(writer, sheet_name="Light_chain", index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-in", "--input_file", help="输入文件")
    parser.add_argument("-out", "--output_file", help="输出文件")
    args = parser.parse_args()

    if args.input_file and args.output_file:
        extract_sequences(args.input_file)
        annotate_seq(args.input_file, args.output_file)
        # annotate_seq("output_vh.fasta", args.output_file)
    else:
        print("请输入输入文件和输出文件的路径")

if __name__ == "__main__":
    main()
