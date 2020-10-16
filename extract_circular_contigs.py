"""
Author: Timo Lucas
This python scripts takes as input a .gfa file from an assembly as well as a file with contigs that were corrected after assembly
it then searches for circular contigs in the .gfa file and extracts the corresponding corrected contigs from the provided assembly
the script also removes the contigs from the original bins in order to prevent redundancy
This script was written to improve the reference binning of metagenomic assemblies by creating new bins for ciruclar contigs
from bins that contain many ambigous contigs

"""

import argparse
import traceback

import Bio.SeqIO as bioseq
import gfapy as gfa

parser = argparse.ArgumentParser(description='Extract circular contigs from metagenome binning, creates new bin with pattern [old name + circular] and removes them from old bin')
# reads_to_reference file created from remap_reads_assembly.sh
parser.add_argument('-g', type=str,
                    help='Path to assembly graph in gfa format', dest='g')
parser.add_argument('-a', type=str,
                    help='Path to corrected assembly in fasta format e.g. assembly after frame shift correction using nanopore assembly pipeline', dest='a')

parser.add_argument('-n', type=str,
                    help='Path to file containing contig annotations (read to class exracted from MEGAN using daa2info tool)')

parser.add_argument('-o', type=str,
                    help='Path to folder with assembly bins. The script saves the new bins here and deletes the contigs from old bins found here.')


args = parser.parse_args()




gfa_path = args.g
assembly_path = args.a
name_file_path = args.n
bin_folder = args.o

circular = []
circular_filtered = []
segment_len_dict = {}
name_dict = {}


contig_names = open(name_file_path,'r')



#load assembly graph
assembly_graph = gfa.Gfa.from_file(gfa_path)

# get name of edges that are connected to themselves


for edge in assembly_graph.edges:
    node1 = str(edge).split('\t')[1]
    node2 = str(edge).split('\t')[3]
    orientation1 = str(edge).split('\t')[2]
    orientation2  = str(edge).split('\t')[4]
    if node1 == node2 and orientation1 == orientation2:
        circular.append(node1)

# get all segment sequences in a dictionary

for segment in assembly_graph.segments:
    name = str(segment).split('\t')[1]
    seq = str(segment).split('\t')[2]
    segment_len_dict[name] = len(seq)

#for all names in circular list remove the ones with len < 200kbp (possible plasmids)


for contig in circular:
    circular_filtered.append(contig)


for i, contig in enumerate(circular_filtered):
    if contig.startswith('edge'):
        circular_filtered[i] = contig.replace('edge','contig')


print(circular_filtered)


#load contig annotation for bin naming
for line in contig_names:
    name = line.split('\t')[0].split('_segment')[0]
    print(name)
    if name in circular_filtered:
        na = str(line.split('\t')[1].replace(' ','_'))
        #QUICK FIX FOR ERROR WITH CAG:878 (there is : after CAG so I replaced it to fit the bin names
        if "CAG" in na:
            na = str(na.replace(':', '_'))
        if "NBRC" in na:
            na = str(na.replace('=', '_'))
        name_dict[name] = na.rstrip("\n")


print(circular_filtered)










#now input corrected assembly contigs
assembly_corrected = bioseq.parse(assembly_path,'fasta')

#grab the circular contigs from the corrected assembly and export it to new bins
for record in assembly_corrected:
    name = record.name.split('_segment')[0]
    if name in circular_filtered:
        try:
            print(f"{bin_folder}{name_dict[name]}.fasta")
            old_bin = bioseq.parse(f"{bin_folder}/{name_dict[name]}.fasta", "fasta")
            old_bin_updated = []

            for old_record in old_bin:
                if old_record.name.split('_segment')[0] not in circular_filtered:
                    old_bin_updated.append(old_record)
            #write updated old bin without the circular contig
            print(f'Deleting contig {name} from old bin {name_dict[name]}.fasta')
            bioseq.write(old_bin_updated,f"{bin_folder}/{name_dict[name]}.fasta",'fasta')
            #write new circular bin
            print(f'Exporting circular contig {name} to new bin called {name_dict[name]}_circular.fasta')
            bioseq.write(record,f'{bin_folder}/{name_dict[name]}_circular.fasta','fasta')
        except:
            traceback.print_exc()
            print(f"No annotation found for contig {name} writing to unknown.fasta")

            try:
                unknown_old = []
                unknown = bioseq.parse(f'{bin_folder}/unknown.fasta', 'fasta')
                for record_u in unknown:
                    unknown_old.append(record_u)
                unknown_old.append(record)
            except:
                print("Unknown empty")
            bioseq.write(unknown_old, f'{bin_folder}/unknown.fasta', 'fasta')



