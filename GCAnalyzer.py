from statistics import median
from Bio import SeqIO
import sys
import csv
from typing import Union


OUT_PATH = "O1P_SeqNumberTable.txt"
OUT_PATH_2 = "O2P_GCcontent.txt"

"""
O1P_SeqNumberTable.txt - this file should contain a tab separated table of the sequence number of Flowers, Animals, and Bacteria
"""

"""
O2P_GCcontent.txt - this file should contain information about the overall GC-content (considering all sequences), and the range and median 
"""



def write_csv(dest,data):
    fd = open(dest, mode="w")
    csv_writer = csv.writer(fd, delimiter=f"\t", quotechar="|")

    for row in data:
        csv_writer.writerow(row)


def read_file(file_path:str):
    fd = open(file_path,"r")
    return SeqIO.parse(fd,"fasta")



#driver function for get_cv_content
def get_gc_info(records:list):
    data = []
    for r in records:
        data.append(get_cv_content(str(r.seq))) #r.seq is the full sequence not truncated
    
    #combine all the sequencies into one
    full_seq = [str(r.seq) for r in records] 
    full_seq = "".join(full_seq)
    
    return round(get_cv_content(full_seq), 3), round(min(data),3), round(max(data), 3), round(median(data),3)


#calculate gc content of a single sequence string
def get_cv_content(seq:str) -> int:
    count_C = seq.count("C")
    count_A = seq.count("A")
    count_G = seq.count("G")
    count_T = seq.count("T")
    
    return (count_G + count_C) / (count_G +count_C + count_T + count_A)



#return amount of a specified seq-name, and a list with that seq  
def count_by_name(name:str,records:list) -> Union[int, list]:
    c = 0
    l = []
    for r in records:
        if name in r.name:
            c += 1
            l.append(r)
    return c,l


def main():
    
    
    if len(sys.argv) < 2:
        print("mangler filsti")
        sys.exit(1)
    
    FILE_PATH = sys.argv[1]
    records = list(read_file(FILE_PATH))

    f_c,f_list = count_by_name("flower", records)
    a_c,a_list = count_by_name("animal",records)
    b_c,b_list = count_by_name("bacterium",records)
    write_csv(OUT_PATH, [["Flowers","Animals","Bacteria"], [f_c,a_c,b_c] ])
    
    f_gc = get_gc_info(f_list)
    b_gc = get_gc_info(b_list)
    a_gc = get_gc_info(a_list)

    write_csv(OUT_PATH_2, 
            [ ["TotalGC FLOWERS:", f_gc[0]], ["MinMaxMedian_Flowers:", f_gc[1], f_gc[2], f_gc[3] ], 
              [],
              ["TotalGC BACTERIA:", b_gc[0]], ["MinMaxMedian_Bacteria:", b_gc[1], b_gc[2], b_gc[3] ],
              [],
              ["TotalGC ANIMALS:", a_gc[0]], ["MinMaxMedian_Animals:", a_gc[1], a_gc[2], a_gc[3] ],
              
              ])


if __name__ == "__main__":
    main()
