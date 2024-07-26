# command: python SCJ_SCF.py <tar> <pseudo_ref> <output_path>

import os
# os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs"
import sys


##############################
# Starting SCJ_SCF.py project
##############################
class HEAD_TAIL_OF_CONTIG_IN_TARGET:
    def __init__(self):
        self.contig_th_map = {}   
        self.th_contig_map = {}  
    def __call__(self, tar_path):
        self.contig_th_map.clear()  # 清空之前的記錄，以確保每次調用都是新的開始
        self.th_contig_map.clear() 
        prev_ctg = ''
        prev_marker = ''
        with open(tar_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                data = line.strip().split()
                ctg_name = data[2]
                marker = data[1]
                if ctg_name != prev_ctg:
                    if prev_ctg != '':
                        self.contig_th_map[prev_ctg].append(prev_marker)  # 幫上一個 ctg 紀錄 tail
                    self.contig_th_map[ctg_name] = [marker]  # 自己紀錄 head
                    prev_ctg = ctg_name
                prev_marker = marker
                if line == lines[-1]: # 最後一行要紀錄最後一個 tail
                    self.contig_th_map[prev_ctg].append(prev_marker)  # 幫上一個 ctg 紀錄 tail

        for ctg, th in self.contig_th_map.items():
            tail = th[0]
            head = th[1]
            self.th_contig_map[tail] = ctg
            self.th_contig_map[head] = ctg
        return self.contig_th_map, self.th_contig_map
    
    # For testing
    def print(self):
        for ctg, ht in self.contig_th_map.items():
            print(f'CTG: {ctg}, TAIL: {ht[0]}, HEAD: {ht[1]}')
        print(f'#TOTAL_CTGS: {len(self.contig_th_map)}')

class BEST_SCJ_ADJ:
    def __init__(self):
        self.best_scj_adj = set()
    def __call__(self, ref_path):
        self.best_scj_adj.clear()  # 清空之前的記錄，以確保每次調用都是新的開始
        prev_ctg = ''
        prev_marker = ''
        with open(ref_path, 'r') as file:
            lines = file.readlines()
            for line in lines:
                data = line.strip().split()[1:]
                for idx in range(len(data)):
                    if idx == len(data)-1:
                        break
                    self.best_scj_adj.add(f"{data[idx]} {data[idx+1]}")
        return self.best_scj_adj
    # For testing
    def print(self):
        print('BEST_SCJ_ADJ:')
        for adj in self.best_scj_adj:
            print(f'{adj}')
        print(f'#TOTAL_BEST_SCJ_ADJ: {len(self.best_scj_adj)}')

def SCJ_core(CONTIG_TH_MAP, TH_CONTIG_MAP, BEST_SCJ_ADJ):
    contig_th_map = CONTIG_TH_MAP.copy()
    th_contig_map = TH_CONTIG_MAP.copy()
    contig_orien_map = {key: 1 for key in contig_th_map.keys()}
    all_contigs = set(contig_th_map.keys())
    print(f'ALL CONTIGS ARE: {all_contigs}')
    finished_contigs = set()

    all_best_adj = BEST_SCJ_ADJ.copy()


    while(len(all_best_adj) != 0):
        
        ADJ = all_best_adj.pop()
        extrem_1 = ADJ.split()[0]
        extrem_2 = ADJ.split()[1]

        print(f"Now searching for {extrem_1} {extrem_2}")

        if (extrem_1 in th_contig_map) and (extrem_2 in th_contig_map):
            ctg1 = th_contig_map[extrem_1]
            ctg2 = th_contig_map[extrem_2]
        elif (extrem_1 in th_contig_map) and (str(-int(extrem_2)) in th_contig_map):
            extrem_2 = str(-int(extrem_2))
            ctg1 = th_contig_map[extrem_1]
            ctg2 = th_contig_map[extrem_2]
        elif (str(-int(extrem_1)) in th_contig_map) and (extrem_2 in th_contig_map):
            extrem_1 = str(-int(extrem_1))
            ctg1 = th_contig_map[extrem_1]
            ctg2 = th_contig_map[extrem_2]
        elif (str(-int(extrem_1)) in th_contig_map) and (str(-int(extrem_2)) in th_contig_map):
            extrem_1 = str(-int(extrem_1))
            extrem_2 = str(-int(extrem_2))
            ctg1 = th_contig_map[extrem_1]
            ctg2 = th_contig_map[extrem_2]
        else:
            print(f"Not find {extrem_1} {extrem_2}!")
            continue
        

        # 若 ctg1 和 ctg2 是同一個，則不進行下面四種可能的 scaffold 動作
        if(ctg1 == ctg2):
            continue

        if (f"{-int(contig_th_map[ctg1][0])} {contig_th_map[ctg2][0]}" in BEST_SCJ_ADJ) or (f"{-int(contig_th_map[ctg2][0])} {contig_th_map[ctg1][0]}" in BEST_SCJ_ADJ):
            for ctg in ctg1.split():
                contig_orien_map[ctg] = -contig_orien_map[ctg]
            rev_ctg1  = ' '.join(reversed(ctg1.split()))

            contig_th_map[f"{rev_ctg1} {ctg2}"] = [str(-int(contig_th_map[ctg1][1])), contig_th_map[ctg2][1]]

            del th_contig_map[extrem_1]
            del th_contig_map[extrem_2]
            if contig_th_map[ctg1][1] in th_contig_map:
                del th_contig_map[contig_th_map[ctg1][1]]
            th_contig_map[str(-int(contig_th_map[ctg1][1]))] = f"{rev_ctg1} {ctg2}"
            th_contig_map[contig_th_map[ctg2][1]] = f"{rev_ctg1} {ctg2}"

            del contig_th_map[ctg1]
            del contig_th_map[ctg2]
            

        
        elif (f"{-int(contig_th_map[ctg1][0])} {-int(contig_th_map[ctg2][1])}" in BEST_SCJ_ADJ) or (f"{contig_th_map[ctg2][1]} {contig_th_map[ctg1][0]}" in BEST_SCJ_ADJ):
            contig_th_map[f"{ctg2} {ctg1}"] = [contig_th_map[ctg2][0], contig_th_map[ctg1][1]]

            del th_contig_map[extrem_1]
            del th_contig_map[extrem_2]
            th_contig_map[contig_th_map[ctg2][0]] = f"{ctg2} {ctg1}"
            th_contig_map[contig_th_map[ctg1][1]] = f"{ctg2} {ctg1}"

            del contig_th_map[ctg1]
            del contig_th_map[ctg2] 
            

        
        elif (f"{contig_th_map[ctg1][1]} {contig_th_map[ctg2][0]}" in BEST_SCJ_ADJ) or (f"{-int(contig_th_map[ctg2][0])} {-int(contig_th_map[ctg1][1])}" in BEST_SCJ_ADJ):
            contig_th_map[f"{ctg1} {ctg2}"] = [contig_th_map[ctg1][0], contig_th_map[ctg2][1]]
            
            del th_contig_map[extrem_1]
            del th_contig_map[extrem_2]  
            th_contig_map[contig_th_map[ctg1][0]] = f"{ctg1} {ctg2}"
            th_contig_map[contig_th_map[ctg2][1]] = f"{ctg1} {ctg2}"

            del contig_th_map[ctg1]
            del contig_th_map[ctg2] 
            

        
        elif (f"{contig_th_map[ctg1][1]} {-int(contig_th_map[ctg2][1])}" in BEST_SCJ_ADJ) or (f"{contig_th_map[ctg2][1]} {-int(contig_th_map[ctg1][1])}" in BEST_SCJ_ADJ):
            for ctg in ctg2.split():
                contig_orien_map[ctg] = -contig_orien_map[ctg]
            rev_ctg2 = ' '.join(reversed(ctg2.split()))                
            
            contig_th_map[f"{ctg1} {rev_ctg2}"] = [contig_th_map[ctg1][0], str(-int(contig_th_map[ctg2][0]))]

            del th_contig_map[extrem_1]
            del th_contig_map[extrem_2]
            if contig_th_map[ctg2][0] in th_contig_map:
                del th_contig_map[contig_th_map[ctg2][0]]
            th_contig_map[contig_th_map[ctg1][0]] = f"{ctg1} {rev_ctg2}"
            th_contig_map[str(-int(contig_th_map[ctg2][0]))] = f"{ctg1} {rev_ctg2}"

            del contig_th_map[ctg1]
            del contig_th_map[ctg2]
            

        else:
            print(f"{ctg1} CAN'T MATCHED WITH {ctg2}!")
    
    finished_contigs = set(contig_th_map.keys())
    
    return finished_contigs, contig_orien_map
            
                

# Main function
def main():

    # python FNA2ALLm.py <# of ref genomes> <tar> <ref1> <ref2> ... <output_path>
    # <output_path> should include the outputs of Sibelia.

    try:
        if len(sys.argv) != 4:
            raise ValueError("Invalid number of arguments")
    except ValueError as e:
        print(f"Usage: python SCJ_SCF.py <tar> <pseudo_ref> <output_path>\nError: {e}")
        sys.exit()

    tar_path = sys.argv[1]
    ref_path = sys.argv[2]
    out_path = sys.argv[3]

    # 爬取 target 中的 ctg 的頭跟尾
    preprocess_1 = HEAD_TAIL_OF_CONTIG_IN_TARGET()
    contig_th_map, th_contig_map = preprocess_1(tar_path)
    preprocess_1.print()

    # 爬取 Best_SCJ_Adjacency 中的最佳相鄰
    preprocess_2 = BEST_SCJ_ADJ()
    best_scj_adj = preprocess_2(ref_path)
    preprocess_2.print()

    # 做 SCJ-scaffolding main part
    finished_contigs, contig_orien_map = SCJ_core(contig_th_map, th_contig_map, best_scj_adj)
    print("FINISH SCAFFOLDING !")
    print("FINAL CONTIGS:")
    for ctg in finished_contigs:
        print(ctg)
    # # 輸出成 ScaffoldResult
    with open(out_path+'/ScaffoldResult', 'w') as file:
        ctg_idx = 0
        for scf in finished_contigs:
            ctgs = scf.split()
            if len(ctgs) == 1:
                continue
            ctg_idx = ctg_idx + 1            
            file.write(f'>scaffold_{ctg_idx}\n')
            for ctg in ctgs:
                orient = int(contig_orien_map[ctg]==-1)
                file.write(f'{ctg} {orient}\n')
            file.write('\n')
            
    
    print("SCJ_SCF.py done!")
    

if __name__ == '__main__':
    main()