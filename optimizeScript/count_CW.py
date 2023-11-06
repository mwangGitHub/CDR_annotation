'''
This script is to count the number of C or W in 1031_ab_sequences (provided by **)
and also calculated the distance beteen C1, W and C2.

In IMGT annotation method, 
	there should have 2 Cys (in position 23-first and 104-second) 
	and 1 W (in position 41) at least in ab sequence.

Results:
in total 120 ab sequences: there are
 26 sequences not qualify the C/W number
 75 sequences not qualify the distance range
 19 sequences qualify the C/W number and distance requirement

while, only 2 of 19 can be annotated by ANARCI

Here are the distance range requirement from C1, W, C2
    
distance_range = {
    'c1c2': (60, 82),
    'c1w':  (10, 19),
    'wc2':  (50, 64)
    }


'''

# 统计C/W的个数及位置信息，输出不满足个数要求的序列及其信息，返回满足个数要求的序列及其信息
def get_CW_info(seq, bad_num=0): 
    # idx_sequence += 1
    C_positions = [i+1 for i in range(len(seq)) if seq[i]=='C']
    W_positions = [i+1 for i in range(len(seq)) if seq[i]=='W']
    num_C = len(C_positions)
    num_W = len(W_positions)
    if num_C < 2 or num_W < 1:
        bad_num += 1
        print(f'the number of C/W is not qualified')
        # print(f'In Sequence: {seq}\t length: {len(seq)}')
        print(f'number of Cys: {len(C_positions)}\npositions: {C_positions}')
        print(f'number of W: {len(W_positions)}\npositions: {W_positions}\n')
        return None, None
    else: 
        # print(f'In Sequence: {seq}\t length: {len(seq)}')
        return C_positions, W_positions
    
def cal_distance(C_list, W_list): # 参数是位置列表
    C1C2_distance = [C_list[i] - C_list[i-1] for i in range(1, len(C_list))]
    C1W_distance = []
    WC2_distance = []

    for i in range(len(C_list)-1):
        for j in range(len(W_list)):
            c1w = W_list[j] - C_list[i]
            wc2 = C_list[i+1] - W_list[j]
            C1W_distance.append(c1w)
            WC2_distance.append(wc2)
    return C1C2_distance, C1W_distance, WC2_distance
        

def in_range(c1c2, distance_range):
    for dis in c1c2:
        if distance_range[0] <= dis <= distance_range[1]:
            return True
    print(f'The distance between {c1c2} not qualify distance range {distance_range}\n')
    return False


def main(input_fasta_file): # 输入fasta文件

    distance_range = {
        'c1c2': (60, 82),
        'c1w':  (10, 19),
        'wc2':  (50, 64)
    }

    with open(input_fasta_file, 'r') as in_file:
        for line in in_file:
            if  line.startswith('>Sequence'):
                pass
                # print(f'In {line}')
            else:
                print(f'*****Sequence: {line}******')
                C_positions, W_positions = get_CW_info(line)
                if  C_positions == None and  W_positions == None:
                    pass
                else:
                    c1c2_dis, c1w_dis, wc2_dis = cal_distance(C_positions, W_positions)
                    if in_range(c1c2_dis, distance_range['c1c2']) and in_range(c1w_dis, distance_range['c1w']) and in_range(wc2_dis, distance_range['wc2']):
                        print(f'qualify the number of C/W and the distance between them')
                        print(f'C positions: {C_positions}\nW positions: {W_positions}\n')
                        print(f'c1c2_dis: {c1c2_dis}\nc1w_dis: {c1w_dis}\nwc2_dis: {wc2_dis}\n')
                    else:
                        print(f'C_positions: {C_positions}\nW_positions: {W_positions}')
                        print(f'c1c2_distance: {c1c2_dis}\n c1w_distance: {c1w_dis} \n wc2_distance: {wc2_dis}\n\n')


if __name__ == '__main__':

    input_fasta_file = '/Users/yscao/mmeng/deeplearning/ANARCI/CDR_annotation/optimizeScript/fasta_ab_1031.fasta'

    main(input_fasta_file)
