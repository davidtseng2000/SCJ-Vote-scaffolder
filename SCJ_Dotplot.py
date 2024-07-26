import matplotlib.pyplot as plt
import argparse


contig_to_marker = {}
tar_ctgs_separate = []
ref_ctgs_separate = []


def read_file(filename, tar_or_ref):
    global tar_ctgs_separate
    if tar_or_ref == 'tar':
        tar_ctgs_separate = []
    data = []
    prev_ctg = ""
    order = 1
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            order = int(parts[0])
            marker = int(parts[1])
            contig = parts[2]
            data.append((order, marker, contig))

            if tar_or_ref == 'tar':
                if contig not in contig_to_marker:
                    contig_to_marker[contig] = [marker]
                    tar_ctgs_separate.append(order)
                else:
                    contig_to_marker[contig].append(marker)
                    order += 1
            
            if tar_or_ref == 'ref':
                if contig != prev_ctg:
                    ref_ctgs_separate.append(order)
                    prev_ctg = contig

    return data

def read_scf(filename):
    global contig_to_marker
    global tar_ctgs_separate
    tar_ctgs_separate = []
    data = []
    order = 1
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) != 2:
                if len(parts) >= 1 and parts[0].startswith('>'):
                    tar_ctgs_separate.append(order)
                continue

            contig = parts[0]

            for marker in contig_to_marker.get(contig, []):
                data.append((order, marker, contig))
                order += 1
    return data


def draw(tar_data, ref_data, outputfile, scf_or_not):

    global tar_ctgs_separate
        
    plt.figure(figsize=(10, 10))


    for t_order, t_marker, t_contig in tar_data:
        for r_order, r_marker, r_contig in ref_data:
            if abs(t_marker) == abs(r_marker):
                plt.plot(r_order, t_order, 'bo') 


    for i, tar_scf in enumerate(tar_ctgs_separate, start=1):
        plt.axhline(y=tar_scf, color='r', linestyle='--')
    for i, ref_ctg in enumerate(ref_ctgs_separate, start=1):
        plt.axvline(x=ref_ctg, color='r', linestyle='--')
        

    # 設定軸標籤
    plt.xlabel('Ref Marker Order')
    plt.ylabel('Target Marker Order')

    if scf_or_not == 0:
        plt.xticks(ticks=ref_ctgs_separate, labels=[f'ctg_{i+1}' for i in range(len(ref_ctgs_separate))])
        plt.yticks(ticks=tar_ctgs_separate, labels=[f'ctg_{i+1}' for i in range(len(tar_ctgs_separate))])
        plt.setp(plt.gca().get_xticklabels(), rotation=90, ha='right')

    if scf_or_not == 1:
        plt.xticks(ticks=ref_ctgs_separate, labels=[f'ctg_{i+1}' for i in range(len(ref_ctgs_separate))])
        plt.yticks(ticks=tar_ctgs_separate, labels=[f'scf_{i+1}' for i in range(len(tar_ctgs_separate))])
        plt.setp(plt.gca().get_xticklabels(), rotation=90, ha='right')

    # 顯示圖表
    plt.title('Marker Match Plot')
    plt.grid(True)

    plt.savefig(outputfile)

    # # 顯示圖表
    # plt.show()

def main():
    parser = argparse.ArgumentParser(description="dotplot program for simulated data")
    parser.add_argument("outputfile", help="outputfile")
    args = parser.parse_args()

    # 讀取 target.all 和 ref.all 檔案
    tar_data = read_file(args.outputfile + '/tar.nondup.all', 'tar')
    ref_data = read_file(args.outputfile + '/Best_SCJ.all', 'ref')
    draw(tar_data, ref_data, args.outputfile+'/DotplotBefore.png', 0)

    scf_data = read_scf(args.outputfile + '/ScaffoldResult')        
    draw(scf_data, ref_data, args.outputfile+'/DotplotAfterSCF.png', 1)

if __name__ == "__main__":
    main()
