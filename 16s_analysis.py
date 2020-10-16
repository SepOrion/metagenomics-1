#script gets as input a folder containing 16s sequences that were aligned against NCBI-nr 16s database as well as
#reads aligned to those sequences. script creates summary file showing 16s sequences, their taxonomic assignment from blast
# and the amount of reads that aligned to those sequences

import glob
import traceback

from Bio import SeqIO

#dictionary (file_name: (contig_name:contig_properties))
rna_dict = {}
contig_dict = {}
out_path = "/Users/timolucas/Documents/jeon/r1_august/16S_summary.tsv"
out_folder = "/Users/timolucas/Documents/jeon/r1_august/"
rna_sequences_dict = {}
contig_to_bin_name_dict = {}
output_file = open(out_path,'w')
tab='\t'
newline='\n'
output_file.write(f"16S RNA sequence {tab}"
                  f" length {tab}"
          f" contig origin {tab}"
          f" number of reads aligned {tab}"
          f" taxonomic assignment by assembly bin {tab}"
          f" taxonomic assignment by 16S{newline}")



output_file.close()

def write_summary(bin_name):
    rna_count_dict = {}


    output_file = open(out_path,'a')
    rna_file = bin_name
    rna_top_hits = open(rna_file + "_top_hits", "r")

    for line in rna_top_hits:
        if line.split("\t")[0].startswith("16S"):
            rna_name = line.split('\t')[0]
            percent_identity = (line.split('\t')[2])
            rna_top_hit= line.split('\t')[8]
            contig_dict[rna_name] = (rna_top_hit,percent_identity)

    rna_read_count_file = open(rna_file+"_reads_aligned",'r')
    read_count_list = []
    for line in rna_read_count_file:
        read_count_dict = {}

        if line.split('\t')[5].startswith('16S'):
            read_count_list.append(line.split('\t')[5])


    rna_read_id_list = list(set(read_count_list))

    for rna in list(set(read_count_list)):
        read_count_dict[rna] = read_count_list.count(rna)



    rna_sequences_file = open(rna_file+'.fasta','r')


    count=0
    for record in SeqIO.parse(rna_file+".fasta", "fasta"):
        if record.id.startswith("16S"):
            count+=1
            rna_sequences_dict[record.id] = record.seq
        rna_count_dict[record.id] = count


#create tab separated output and add all important information to it

    for rna in rna_read_id_list:
        contig_to_bin_name_dict[rna] = bin_name.split('/')[-1]
        output_file.write(f"{rna_sequences_dict[rna]} {tab}"
              f" {len(rna_sequences_dict[rna])} {tab}"
              f" {rna} {tab}"
              f" {read_count_dict[rna]} and {read_count_dict[rna]/rna_count_dict[rna] }  {tab}"
              f" {rna_file.split('/')[-1]} {tab}"
              f" {contig_dict[rna]}{newline}")

#also create one file containing all 16s sequences from the sample







for contig_bin in glob.glob("/Users/timolucas/Documents/jeon/r1_august/extracted_reads/*_rna.fasta"):
    try:
        write_summary(contig_bin.split('.fasta')[0])
    except:
         print(f"Exception occured while writing {contig_bin}")
         traceback.print_exc()


out = open(f"{out_folder}/16S_sequences.fasta",'w')

print(contig_to_bin_name_dict)
try:
    for key in rna_sequences_dict.keys():
        out.write(f">{contig_to_bin_name_dict[key]}\n")
        out.write(f"{rna_sequences_dict[key]}\n")
except:
    print(f"exception occurred with while writing {key}")
    traceback.print_exc()
