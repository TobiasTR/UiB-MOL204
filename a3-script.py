
from Bio import SeqIO
import argparse
import os


def read_file(file_path:str):
    fd = open(file_path,"r")
    l =  list(SeqIO.parse(fd,"fasta"))
    return l


def write_fasta(data:list, path:str):
    fd = open(path,"w")
    for r in data:
        SeqIO.write(r,fd,"fasta")
    fd.close()


def sort_by_len(records):
    return sorted(records, key= lambda x: len(x.seq), reverse=True)


def dna_to_rna(seq):
    seq.seq = seq.seq.complement_rna()
    return seq


def reverse_complement(seq):
    seq.seq = seq.seq.reverse_complement()
    return seq

#transcribe the dna to rna then translate to a protein sequence
def translate(seq):
    seq.seq = seq.seq.transcribe()
    seq.seq = seq.seq.translate(table=1)
    return seq


def main():

    parser = argparse.ArgumentParser(description='read fasta file, create new file with new output in current directory.')
    parser.add_argument('-i','--input', help='path to fasta file',required=True)
    parser.add_argument('-dtr','--dna2rna', help='convert to rna',action='store_true')
    parser.add_argument('-rc','--reverse', help='create reverse complement',action='store_true' )
    parser.add_argument('-tr','--translate', help='translate mRNA to amino acid sequence',action='store_true' )
    parser.add_argument('-so','--sort', help='sort by seq len',action='store_true' )
    args = parser.parse_args()    
    dir_path = os.path.dirname(os.path.realpath(__file__)) + "/"

    FILE_PATH = args.input
    records = read_file(FILE_PATH)

    if args.dna2rna:
        r = list(map(dna_to_rna, records))
        write_fasta(r, dir_path + "rna.fasta")

    if args.reverse:
        r = list(map(reverse_complement, records))
        write_fasta(r,dir_path + "reverse.fasta")
    
    if args.translate:
        r = list(map(translate, records))
        write_fasta(r,dir_path + "proteins.fasta")
    
    if args.sort:
        write_fasta(sort_by_len(records), dir_path + "sorted.fasta")

    
    return


if __name__ == "__main__":
    main()
