# parse arguments
import argparse

import pandas as pd
import plotly.express as px

parser = argparse.ArgumentParser(description='Plots bacterial abundances')
# reads_to_reference file created from remap_reads_assembly.sh
parser.add_argument('-r2c', type=str,
                    help='2 column file with read name in first column and aligned contig in 2nd col', dest='r2c')
parser.add_argument('-c2t', type=str,
                    help='2 column file with contig name and assigned MEGAN taxonomy in 2nd column', dest='c2t')

parser.add_argument('-reads', type=str,
                    help='one column file with read names (I use this to search for the reads that did not align to pseudomonas)')

args = parser.parse_args()


# load read to contig and contig to taxonomy file


# reads_file = open(args.reads,'r')

def create_plot(r2c, c2t, plot_title, plot_name):
    print(r2c)
    print(c2t)
    taxonomy_dict = {}
    annotated_contigs = []
    contig_count_dict = {}
    read_to_contig_dict = {}
    read_to_taxonomy_dict = {}

    r2c = open(r2c, 'r')
    c2t = open(c2t, 'r')
    for line in c2t.read().splitlines():
        contig_name = line.split("\t")[0]
        contig_taxonomy = line.split("\t")[1]
        taxonomy_dict[contig_name] = contig_taxonomy
    count = 0
    for line in r2c.read().splitlines():
        try:
            read = line.split("\t")[0]
            aligned_contig = line.split("\t")[1]
            annotated_contigs.append(taxonomy_dict[aligned_contig])
            read_to_taxonomy_dict[read] = taxonomy_dict[aligned_contig]
        except KeyError:
            count += 1

    # pseudomonas_reads_better=[]
    # for line in reads_file.read().splitlines():
    #     try:
    #         pseudomonas_reads_better.append(read_to_taxonomy_dict[line])
    #     except:
    #         print("key error")
    # print(pseudomonas_reads_better)
    # print(pseudomonas_reads_better.count("unclassified Clostridiales (miscellaneous)"))

    # print(count)
    annotated_contigs_dedup = list(set(annotated_contigs))

    for contig in annotated_contigs_dedup:
        contig_count_dict[contig] = annotated_contigs.count(contig)

    print(contig_count_dict)

    df = pd.DataFrame(contig_count_dict.items(), columns=['species', 'count'])
    print(df)
    fig = px.pie(df, values='count', names='species', title=plot_title)
    # fig.to_image(format="", engine="kaleido")
    # fig.write_image("/Users/timolucas/Documents/jeon/remap_reads/"+plot_name)
    # fig = px.bar(df, x="species", y="count", color='count', barmode='group',
    #             height=400)

    fig.show()


# create_plot(args.r2c, args.c2t, "Abundances R1 August20", "abundances_R1August20.html")
# set plotly orca executable
# plt.io.orca.config .executable = '/Users/timolucas/miniconda3/bin/orca'
# create_plot("/Users/timolucas/Documents/jeon/remap_reads/2n3_reads_to_ref.txt","/Users/timolucas/Documents/jeon/remap_reads/2n3_contigs_to_gtdb.txt","Jeon sample 2n3 abundance plot using gtdb binning","2n3_gtdb.svg")
# create_plot("/Users/timolucas/Documents/jeon/remap_reads/2n3_reads_to_ref.txt","/Users/timolucas/Documents/jeon/remap_reads/2n3_contigs_to_ncbi.txt","Jeon sample 2n3 abundance plot using ncbi binning","2n3_ncbi.png")

# create_plot("/Users/timolucas/Documents/jeon/remap_reads/7n8_reads_to_ref.txt","/Users/timolucas/Documents/jeon/remap_reads/7n8_contigs_to_gtdb.txt","Jeon sample 7n8 abundance plot using gtdb binning","7n8_gtdb.html")
# create_plot("/Users/timolucas/Documents/jeon/remap_reads/7n8_reads_to_ref.txt","/Users/timolucas/Documents/jeon/remap_reads/7n8_contigs_to_ncbi.txt","Jeon sample 7n8 abundance plot using ncbi binning","7n8_ncbi.png")

# create_plot("/Users/timolucas/Documents/jeon/remap_reads/9_reads_to_reference_filtered.txt","/Users/timolucas/Documents/jeon/remap_reads/9_contigs_to_gtdb.txt","Jeon sample 9 abundance plot using gtdb binning","9_gtdb.html")
# create_plot("/Users/timolucas/Documents/jeon/remap_reads/9_reads_to_reference_filtered.txt","/Users/timolucas/Documents/jeon/remap_reads/9_contigs_to_ncbi.txt","Jeon sample 9 abundance plot using ncbi binning","9_ncbi.png")

# create_plot("/Users/timolucas/Documents/jeon/remap_reads/7n8_new_gtdb_reads_to_ref.txt","/Users/timolucas/Documents/jeon/remap_reads/7n8_new_gtdb_contigs_to_taxonomy.txt","Jeon sample 7n8 abundance plot using new gtdb mapping file","7n8_new_gtdb.html")

create_plot('/Users/timolucas/Documents/jeon/r1_august/r1_august_reads_to_reference.txt','/Users/timolucas/Documents/jeon/r1_august/r1_august_contigs_to_taxonomy.txt',"Jeon sample August_R1 abundance plot using new gtdb mapping file","/Users/timolucas/Documents/jeon/r1_august/R1_august.html")
create_plot('/Users/timolucas/Documents/jeon/r2_august/r2_august_reads_to_reference.txt','/Users/timolucas/Documents/jeon/r2_august/r2_august_contig_to_class.txt',"Jeon sample August_R2 abundance plot using new gtdb mapping file","/Users/timolucas/Documents/jeon/r2_august/R2_august.html")

