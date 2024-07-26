# command: python FNA2ALLm.py <# of ref genomes> <tar> <ref1> <ref2> ... <output_path>

import os
# os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs"
import sys
import itertools

# Global variable
ref_num = -1


##############################
# Starting SCJ_Vote_Sim.py project
##############################
genome_contig_map = {}
contig_genome_map = {}
def genome_contigs(file_path, genome_name):
    global genome_contig_map
    global contig_genome_map
    genome_contig_map[genome_name] = []
    prev_ctg = ""
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()
            ctg_name = line.split()[2]  
            if(ctg_name != prev_ctg):
                genome_contig_map[genome_name].append(ctg_name)
                prev_ctg = ctg_name

            if ctg_name not in contig_genome_map:
                contig_genome_map[ctg_name] = genome_name

contig_marker_map = {}
def contig_markers(file_path):
    global contig_marker_map
    cur_ctg = ''
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split()
            ctg_name = line[2]  
            marker = int(line[1])
            if ctg_name != cur_ctg:
                cur_ctg = ctg_name
                contig_marker_map[cur_ctg] = [marker]
            else:
                contig_marker_map[cur_ctg].append(marker)




genome_markers_map = {}
genome_specialty_markers_map = {}
contig_specialty_marker_map = {}
valid_markers = set()
def remove_dup_markers(out_path):
    global genome_contig_map
    global contig_marker_map
    global genome_markers_map
    global genome_specialty_markers_map
    global contig_specialty_marker_map
    global valid_markers
    global ref_num

    marker_num = 0
    remove_ctg = 0
    # 第一部分: 計算各個 genome (tar or ref1 or ref2 or ...) 中的 marker 數，刪掉重複者
    for g, ctgs in genome_contig_map.items():
        # 計算某個 genome (tar or ref1 or ref2 or ...) 中的 marker 數
        marker_cnt = {}
        for ctg in ctgs:
            # Sibelia 的 genomes_permutations 結果 (contig_marker_map) 中可能沒有包含所有 ctgs 
            if ctg not in contig_marker_map:
                print(f'Do NOT have {ctg} data!')
                continue
            for m in contig_marker_map[ctg]:
                if abs(m) not in marker_cnt:
                    marker_cnt[abs(m)] = 1
                else:
                    marker_cnt[abs(m)] = marker_cnt[abs(m)]+1 
        # 刪除重複出現的 marker
        non_dup_markers = set()
        for m, cnt in marker_cnt.items():
            if cnt == 1:
                non_dup_markers.add(m)
        # 更新 contig_marker_map，僅留下沒重複出現的 marker
        for ctg in ctgs:
            if ctg not in contig_marker_map:
                continue
            contig_marker_map[ctg] = [m for m in contig_marker_map[ctg] if abs(m) in non_dup_markers]
            if len(contig_marker_map[ctg]) == 0:
                remove_ctg = remove_ctg + 1
            else:
                marker_num = marker_num + len(contig_marker_map[ctg])
        genome_markers_map[g] = non_dup_markers
    with open(str(out_path)+'/contig_marker_analysis.txt', 'a') as file:           
        file.write(f'\n')
        file.write(f'After remove_part1:\n')
        file.write(f'Remove #ctg: {remove_ctg}\n')
        file.write(f'#All markers: {marker_num}\n')

    # 第二部分: 且只留下在各個 genome (tar or ref1 or ref2 or ...) 中，恰好出現一次的 marker
    marker_num = 0
    remove_ctg = 0
    all_genomes = ['tar']
    for idx in range(ref_num):
        all_genomes.append(f'ref{idx+1}')
    marker_sets = [genome_markers_map[genome] for genome in all_genomes]
    valid_markers = set.intersection(*marker_sets)
    # valid_markers = (genome_markers_map['tar'] & genome_markers_map['ref1'] & genome_markers_map['ref2'] & genome_markers_map['ref3'])
    print(f'VALID_MARKERS: {valid_markers}')
    # # 240313 新增，另外紀錄 "僅出現在某個 ref 中，而未出現在其他 ref 中的 specialty markers"
    # genome_specialty_markers_map['ref1'] = (genome_markers_map['ref1']-(genome_markers_map['ref2'] | genome_markers_map['ref3']) )
    # genome_specialty_markers_map['ref2'] = (genome_markers_map['ref2']-(genome_markers_map['ref1'] | genome_markers_map['ref3']) )
    # genome_specialty_markers_map['ref3'] = (genome_markers_map['ref3']-(genome_markers_map['ref1'] | genome_markers_map['ref2']) )

    
    for g, ctgs in genome_contig_map.items():
        for ctg in ctgs:
            # Sibelia 的 genomes_permutations 結果 (contig_marker_map) 中可能沒有包含所有 ctgs
            if ctg not in contig_marker_map:
                continue
            # 240313 新增下面兩行
            # if g != 'tar':
            #     contig_specialty_marker_map[ctg] = [m for m in contig_marker_map[ctg] if abs(m) in genome_specialty_markers_map[g]]
            contig_marker_map[ctg] = [m for m in contig_marker_map[ctg] if abs(m) in valid_markers]
            if len(contig_marker_map[ctg]) == 0:
                remove_ctg = remove_ctg + 1
            else:
                marker_num = marker_num + len(contig_marker_map[ctg])
    with open(str(out_path)+'/contig_marker_analysis.txt', 'a') as file:           
        file.write(f'\n')
        file.write(f'After remove_part2:\n')
        file.write(f'Remove #ctg: {remove_ctg}\n')
        file.write(f'#All markers: {marker_num}\n')
    


# Main function
def main():

    # python FNA2ALLm.py <# of ref genomes> <tar> <ref1> <ref2> ... <output_path>
    # <output_path> should include the outputs of Sibelia.

    try:
        if len(sys.argv) != int(sys.argv[1])+4:
            raise ValueError("Invalid number of arguments")
    except ValueError as e:
        print(f"Usage: python SCJ_Vote_Sim.py <# of ref genomes> <tar> <ref1> <ref2> <ref3> ... <output_path>\nError: {e}")
        sys.exit()

    global ref_num
    ref_num = int(sys.argv[1])
    tar_path = sys.argv[2]
    ref_list = sys.argv[3:3+ref_num]
    out_path = sys.argv[-1]
    # Testing
    # print(ref_num)
    # print(tar_path)
    # print(ref_list)
    # print(out_path)
    # print(type(ref_num))
    # print(type(tar_path))
    # print(type(ref_list))
    # print(type(out_path))

    ###############################
    # parse the "genome-contigs" relation
    ###############################
    genome_contigs(tar_path, 'tar')
    for ref_i in range(ref_num):
        ref_name = f'ref{ref_i+1}'
        genome_contigs(ref_list[ref_i], ref_name)
    # Testing
    # for g, ctgs in genome_contig_map.items():
    #     print("-" * 30)
    #     print(g)
    #     for ctg in ctgs:
    #         print(ctg)
        

        
    ###############################
    # parse the "contig-markers" relation
    ###############################
    print('\033[41m' + "=" *50 + '\033[0m') 
    # contig_markers(str(out_path)+'/genomes_permutations.txt')
    contig_markers(tar_path)
    for ref_i in range(ref_num):
        contig_markers(ref_list[ref_i])
    
    # Testing 
    all_ctg_num = 0
    unused_ctg_num = 0
    marker_num = 0   
    print('\033[41m' + 'Part1: Parse the "contig-markers" relation.' + '\033[0m') 
    print('\033[41m' + "=" *50 + '\033[0m') 
    for ctg, markers in contig_marker_map.items():
        try: 
            print("-" * 30)
            print(ctg + ' ' + contig_genome_map[ctg])
            print(markers)
        except:
            print(f"NOBODY KNOWS WHERE {ctg} FROM!")
    # Output
    for g, ctgs in genome_contig_map.items():
        idx = 1
        with open(str(out_path)+'/'+g+'.ori.all', 'w') as file:
            for ctg in ctgs:
                all_ctg_num = all_ctg_num + 1
                if ctg not in contig_marker_map:
                    continue
                if len(contig_marker_map[ctg]) == 0:
                    unused_ctg_num = unused_ctg_num + 1
                for m in contig_marker_map[ctg]:
                    file.write(f'{idx} {m} {ctg} 0\n')
                    idx = idx+1
                    marker_num = marker_num + 1
    with open(str(out_path)+'/contig_marker_analysis.txt', 'w') as file:           
        file.write(f'After parsing:\n')
        file.write(f'#All contigs: {all_ctg_num}\n')
        file.write(f'#Used contigs: {all_ctg_num - unused_ctg_num}\n')
        file.write(f'%Used contigs: {(all_ctg_num - unused_ctg_num)/all_ctg_num}\n')
        file.write(f'#All markers: {marker_num}\n')
            


    ###############################
    # remove duplicated markers in each genome (tar, ref1, ref2, ...)
    ############################### 
    print('\033[41m' + "=" *50 + '\033[0m') 
    print('\033[41m' + 'Part2: Remove duplicated markers in each genome (tar, ref1, ref2, ...).' + '\033[0m')
    print('\033[41m' + "=" *50 + '\033[0m')  
    remove_dup_markers(out_path)   
    # Testing    
    all_ctg_num = 0
    unused_ctg_num = 0
    marker_num = 0
    for ctg, markers in contig_marker_map.items():        
        try:
            print("-" * 30)
            print(ctg + ' ' + contig_genome_map[ctg])
            print(markers)            
        except:
            print(f"NOBODY KNOWS WHERE {ctg} FROM!")
    # Output
    for g, ctgs in genome_contig_map.items():
        idx = 1
        with open(str(out_path)+'/'+g+'.nondup.all', 'w') as file:
            for ctg in ctgs:
                all_ctg_num = all_ctg_num + 1
                if ctg not in contig_marker_map:
                    continue
                if len(contig_marker_map[ctg]) == 0:
                    unused_ctg_num = unused_ctg_num + 1
                for m in contig_marker_map[ctg]:
                    file.write(f'{idx} {m} {ctg} 0\n')
                    idx = idx+1
                    marker_num = marker_num + 1
    with open(str(out_path)+'/contig_marker_analysis.txt', 'a') as file:           
        file.write(f'\n')
        file.write(f'After remove duplications:\n')
        file.write(f'#All contigs: {all_ctg_num}\n')
        file.write(f'#Used contigs: {all_ctg_num - unused_ctg_num}\n')
        file.write(f'%Used contigs: {(all_ctg_num - unused_ctg_num)/all_ctg_num}\n')
        file.write(f'#All markers: {marker_num}\n')


    ###############################
    # calculate the "best" SCJ adjacency
    ############################### 
    # 找各個 ref 中的 adjacency-score
    adj_in_each_ref = []
    for g, ctgs in genome_contig_map.items():
        if 'ref' not in g:
            continue
        adj_score_map = {}        
        for ctg in ctgs:
            if ctg not in contig_marker_map:
                continue
            prev_marker = -1
            for m in contig_marker_map[ctg]:
                if prev_marker == -1:
                    prev_marker = m
                else:
                    if (prev_marker, m) not in adj_score_map:
                        adj_score_map[(prev_marker, m)] = 1
                        adj_score_map[(-m, -prev_marker)] = 1
                    else:
                        adj_score_map[(prev_marker, m)] = adj_score_map[(prev_marker, m)] + 1
                        adj_score_map[(-m, -prev_marker)] = adj_score_map[(-m, -prev_marker)] + 1
                    # print(f'{(prev_marker, m)} in {g}, score:{adj_score_map[(prev_marker, m)]}')
                    prev_marker = m
                    
        adj_in_each_ref.append(adj_score_map)

    # 240717 新增：也去找 target 中的 adjacency
    adj_in_each_tar = []
    for ctg in genome_contig_map['tar']:
        if ctg not in contig_marker_map:
            continue
        prev_marker = -1
        for m in contig_marker_map[ctg]:
            if prev_marker == -1:
                prev_marker = m
            else:
                adj_in_each_tar.append((prev_marker, m))
                prev_marker = m
    with open(str(out_path)+'/adj_in_tar.txt', 'w') as file:           
        for (marker1, marker2) in adj_in_each_tar:
            file.write(f'{marker1} {marker2}\n')                


    # 合併所有 refs 中的 adjacency-score
    adj_score_map_global = {} # 所有 refs 的 adjs
    adj_score_map_best = {} # 投票過半的 adjs
    best_adj_map = {}
    best_adj_set = set()
    best_adj_str_list = []
    for mp in adj_in_each_ref:
        for (marker1, marker2), score in mp.items():
            if (marker1, marker2) not in adj_score_map_global:
                adj_score_map_global[(marker1, marker2)] = score
            else:
                adj_score_map_global[(marker1, marker2)] = adj_score_map_global[(marker1, marker2)] + score
    adj_score_map_global = dict(sorted(adj_score_map_global.items(), key=lambda x:x[1], reverse=True))
    adj_score_map_best = dict(x for x in adj_score_map_global.items() if x[1]>float(ref_num/2))

    with open(str(out_path)+'/all_adj_score.txt', 'w') as file:           
        for (marker1, marker2), score in adj_score_map_global.items():
            file.write(f'{marker1} {marker2} {score}\n')

    
    # Testing 
    print('\033[41m' + "=" *50 + '\033[0m') 
    print('\033[41m' + 'Part3: Making the "adj-score" relation.' + '\033[0m')
    print('\033[41m' + "=" *50 + '\033[0m')
    print("The best adjacency-score obtained from all refs are:")      
    for adj, score in adj_score_map_best.items():
        print("-" * 30)        
        print(f'{adj}, score:{score}')

    # 將 best adjacency 串起來
    for (m1,m2) in adj_score_map_best.keys():
        if m1 not in best_adj_map:
            best_adj_map[m1] = m2            
            best_adj_set.add((m1,m2)) 

    # 建立 (map) best_adj_list
    # 形式如:
    # (key, is a int) -451 -> (value, is a list) [-451, 231, 22, -37, ...]     
    best_adj_lists = {}
    while(len(best_adj_set)!=0):
        m1,m2 = best_adj_set.pop()
        start = m1
        best_adj_lists[start] = [m1, m2]
        best_adj_set.remove((-m2,-m1))
        while m2 in best_adj_map:
            m1 = m2
            m2 = best_adj_map[m1]
            if (m1,m2) in best_adj_set:
                best_adj_lists[start].append(m2)
                best_adj_set.remove((m1,m2))
                best_adj_set.remove((-m2,-m1))
            else:
                # m1 != start to avoid circular
                if m1 in best_adj_lists and (m1 != start):
                    best_adj_lists[start].pop()
                    best_adj_lists[start] = best_adj_lists[start] + best_adj_lists[m1]
                    del best_adj_lists[m1]
                else:                    
                    break
        # dt2000 240624 修正
        # 如果先前已經接了某一段 -650 488 390 -264
        # 後面可能還會接出 650 -249 228 -527 
        # 此時須將後面這段倒過來接到前面那段當中
        if(-start in best_adj_lists):
            best_adj_lists[best_adj_lists[start][-1]] = [-x for x in best_adj_lists[start][::-1]]
            best_adj_lists[best_adj_lists[start][-1]].pop()
            best_adj_lists[best_adj_lists[start][-1]] = best_adj_lists[best_adj_lists[start][-1]] + best_adj_lists[-start]
            del best_adj_lists[start]
            del best_adj_lists[-start]

    # Testing
    # for start, list in best_adj_lists.items():
    #     print('-'*50)
    #     print(start)
    #     print(list)
    #     print('-'*50)
    
    # Output
    with open(str(out_path)+'/'+'Best_SCJ_adjacency.txt', 'w') as file:
        print_idx = 0
        for start_marker, adj_list in best_adj_lists.items():
            print_idx = print_idx + 1
            file.write(f'{print_idx}: ')
            for marker in adj_list:
                file.write(f'{marker} ')
            file.write('\n')


        
    with open(str(out_path)+'/'+'Best_SCJ.all', 'w') as file:
        print_idx = 1
        pseudo_ctg_idx = 0
        for start_marker, adj_list in best_adj_lists.items():
            pseudo_ctg_idx = pseudo_ctg_idx + 1
            for marker in adj_list:
                file.write(f'{print_idx} {marker} ctg_{pseudo_ctg_idx} 1\n')
                print_idx = print_idx + 1
        # 240313 新增
        # for g, ctgs in genome_contig_map.items():
        #     if g == 'tar':
        #         continue
        #     pseudo_ctg_idx = pseudo_ctg_idx + 1
        #     for ctg in ctgs:
        #         # Sibelia 的 genomes_permutations 結果 (contig_marker_map) 中可能沒有包含所有 ctgs
        #         if ctg not in contig_marker_map:
        #             continue
        #         for marker in contig_specialty_marker_map[ctg]:
        #             file.write(f'{print_idx} {marker} ctg_{pseudo_ctg_idx} 1\n')
        #             print_idx = print_idx + 1

    print("FNA2ALLm.py done!")
    

if __name__ == '__main__':
    main()