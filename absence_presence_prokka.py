# This script takes as input a fasta with protein sequences and an alignment (.m8) to create an absence presence plot
# for alignment all bacteria that were confirmed to perform RBOX were downloaded and the enzymes important for RBOX were extracted
# coding sequences for genome assemblies annotated by prokka were then aligned to these proteins to get the alignment file

import glob
import os

import numpy as np
import plotly.figure_factory as ff
from Bio import SeqIO

path_to_rbox_database = "/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa"


def find_proteins_in_gbf_by_ec():
    gbf_file = open("/Users/timolucas/Desktop/rbox_map_test/r1t1_prokka.gbf", "r")
    ec_numbers_rbox = ["1.1.1.1", "1.2.1.10", "2.3.1.9", "1.1.1.35", "1.1.1.157", "1.1.1.35", "1.1.1.157", "1.1.1.35",
                       "4.2.1.17", "7.2.1.2", "3.1.2.21", "3.1.2.20", "4.2.1.55", "2.8.3.8", "2.8.3.9"]

    genes = []
    ec_numbers_gbf = []

    for line in gbf_file:
        if "EC_number" in line:
            ec_number = line.split("=")[1]
            # delete \n
            ec_number = ec_number[:-2]
            ec_number = ec_number[1:]
            ec_numbers_gbf.append(ec_number)

    for ec in ec_numbers_rbox:
        if ec in ec_numbers_gbf:
            print("EC number " + str(ec) + " is in gbf file")


# rbox_identifiers_all_drive = ["K00001 alcohol dehydrogenase [EC:1.1.1.1]", "K04072 acetaldehyde dehydrogenase / alcohol dehydrogenase [EC:1.2.1.10 1.1.1.1]",
#                               "K00626 acetyl-CoA C-acetyltransferase [EC:2.3.1.9]", "K07516 3-hydroxyacyl-CoA dehydrogenase [EC:1.1.1.35]",
#                               "K00074 3-hydroxybutyryl-CoA dehydrogenase [EC:1.1.1.157]",
#                               "K01782 3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase [EC:1.1.1.35 4.2.1.17 5.1.2.3]",
#                               "K01692 enoyl-CoA hydratase [EC:4.2.1.17]",
#                               "K17865 3-hydroxybutyryl-CoA dehydratase [EC:4.2.1.55]",
#                               "K00249 acyl-CoA dehydrogenase [EC:1.3.8.7]", "K00248 butyryl-CoA dehydrogenase [EC:1.3.8.1]",
#                               "K01034 acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]",
#                               "K01034 acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]", "K01073 acyl-CoA hydrolase [EC:3.1.2.20]",
#                               "K01071 medium-chain acyl-[acyl-carrier-protein] hydrolase [EC:3.1.2.21]",
#                               "K00625 phosphate acetyltransferase [EC:2.3.1.8]", "K00925 acetate kinase [EC:2.7.2.1]",
#                               "K03614 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit D [EC:7.2.1.2]",
#                               "K03615 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit C [EC:7.2.1.2]",
#                               "K03616 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit B [EC:7.2.1.2]",
#                               "K00532 ferredoxin hydrogenase [EC:1.12.7.2]","K17998",
#                               "K03522 electron transfer flavoprotein alpha subunit",
#                               "K03521 electron transfer flavoprotein beta subunit",
#                               "K14086 ech hydrogenase subunit A",
#                               "K14087 ech hydrogenase subunit B"]

rbox_protein_names = ["alcohol dehydrogenase",
                      "acetaldehyde dehydrogenase",
                      "acetyl-CoA C-acetyltransferase",
                      "3-hydroxyacyl-CoA dehydrogenase",
                      "3-hydroxybutyryl-CoA dehydrogenase",
                      # "3-hydroxybutyryl-CoA epimerase",
                      "enoyl-CoA hydratase",
                      "3-hydroxybutyryl-CoA dehydratase",
                      "acyl-CoA dehydrogenase",
                      "butyryl-CoA dehydrogenase",
                      "butyryl-CoA:acetate CoA-transferase",
                      "acetyl-CoA hydrolase/transferase",
                      "phosphate acetyltransferase",
                      "acetate kinase",
                      "acetyl-CoA hydrolase",
                      # "acetyl-CoA synthetase",
                      # "trans-2-enoyl-CoA reductase",
                      # "glyceraldehyde-3-phosphate dehydrogenase",
                      # "butyryl-CoA:acetoacetate CoA-transferase"
                      "ferredoxin hydrogenase",
                      "ech hydrogenase"
                      ]

fab_protein_names = ["alcohol dehydrogenase", "aldehyde reductase",
                     "anaerobic carbon-monoxide dehydrogenase", "carbon-monoxide dehydrogenase",

                     "acetyl-CoA decarbonylase",
                     "acetyl-CoA carboxylase",
                     "biotin carboxylase",
                     "S-malonyltransferase",
                     "3-oxoacyl-ACP synthase",
                     "3-oxoacyl-ACP reductase", "3-hydroxyacyl-ACP dehydratase"
    , "enoyl-ACP reductase",
                     "beta-ketoacyl-ACP synthase"]


# clostridiales catherine r1t1
# alignments_to_rbox = open("/Users/timolucas/Documents/spirito/new_annotation/test_rbox_map/r1t1_prokka.m8","r")
# alignments_to_rbox = open("/Users/timolucas/Documents/spirito/assembly/r1t1/extracted_reads/read_id2507162_nameCaproiciproducens_sp._NJN-50/r1t1_prokka.m8","r")


# function gets RBOX fastas database and an alignment to this file as well as proteins that should be considered
# it parses the files and looks the accessions of the alignments up in the original fasta to identify the rbox proteins for that alignments were found
# returns a deduplicated list with names of RBOX proteins
# considered proteins is a filter to limit the search

def get_proteins_in_alignments(fasta_with_rbox_proteins, alignment_file, considered_proteins):
    alignments_to_searched_proteins = open(
        alignment_file,
        "r")

    functional_proteins_aligned = []
    rbox_proteins_aligned_function = []
    # after prokka annotations were aligned vs rbox proteins with following command: time /home3/lucas/software/diamond blastp -d /abprojects/daniel/ncbi/ncbi-Oct2019/nr -q /abprojects/lucas/spirito/new_annotations/rbox_proteins.faa -o /abprojects/lucas/spirito/new_annotations/rbox_proteins_id90_query-cover90.m8 -p 32 -t /dev/shm -c1 -k0 --more-sensitive --taxonmap /abprojects/daniel/ncbi/ncbi-Oct2019/prot.accession2taxid.gz --taxonnodes /abprojects/daniel/ncbi/ncbi-Oct2019/taxdmp.zip --id 80 --query-cover 80 -f 6 qseqid qlen sseqid slen qstart qend sstart send evalue bitscore score length pident positive nident
    # now: go through alignments of prokka CDS to rbox proteins and identify proteins by looking accession up in original RBOX proteins .faa

    # extract all accessions to that alignments were found
    for alignment in alignments_to_searched_proteins:
        # print(alignment)

        line_list = alignment.split("\t")
        protein = line_list[2]
        alignment_length = line_list[11]
        alignment_num_positives = line_list[13]
        percent_positives = int(alignment_num_positives) / int(alignment_length)
        # only consider alignments with at least 50% positive aligned bases
        if percent_positives >= 0.5:
            functional_proteins_aligned.append(protein)

    # print(rbox_proteins_aligned)

    # get the names of the sequences that alignments were found for
    otu_proteins_dict = {}
    for record in SeqIO.parse(fasta_with_rbox_proteins, "fasta"):
        for protein in functional_proteins_aligned:
            if protein == record.name:
                rbox_proteins_aligned_function.append(record.description)
    # deduplicate list
    rbox_proteins_aligned_function = list(dict.fromkeys(rbox_proteins_aligned_function))
    rbox_proteins_aligned_trimmed = []
    # trim accessions

    otu_proteins_dict["Clostridiales"] = rbox_proteins_aligned_function

    # this function gets the .m8 alignment file and parses all

    rbox_proteins_with_alignments = []

    # check if rbox proteins are in clostridiales for alignment file and add them to new list
    # this block of code is not for Clostridiales but for all alignments: TODO: refactor this part of code, works but is confusing
    for protein in considered_proteins:
        for rbox_protein_aligned in otu_proteins_dict["Clostridiales"]:
            if protein in rbox_protein_aligned:
                rbox_proteins_with_alignments.append(protein)
    # deduplicate list with rbox proteins
    rbox_proteins_with_alignments = list(dict.fromkeys(rbox_proteins_with_alignments))
    return (rbox_proteins_with_alignments)


# takes list with proteins computed from get_rbox_proteins_in_alignments and creates heatmap out of it


def create_list_for_plot(proteins_found, all_rbox_proteins):
    # create two lists one for numbers 0 1 and coloring the other for the symbols that are plotted in the end
    data_for_heatmap_one_file = []
    annotation_for_heatmap = []
    for searched_protein in all_rbox_proteins:
        if searched_protein in proteins_found:
            data_for_heatmap_one_file.append(1)
            annotation_for_heatmap.append("X")
        else:
            data_for_heatmap_one_file.append(0)
            annotation_for_heatmap.append(("-"))
    return (data_for_heatmap_one_file, annotation_for_heatmap)


def create_heatmap(plot_list_all_files, bins, heatmap_annotations, protein_names):
    #
    plot_list_all_files = np.array(plot_list_all_files).transpose()
    heatmap_annotations = np.array(heatmap_annotations).transpose()
    # plot_list_all_files = pd.DataFrame(plot_list_all_files)
    print(plot_list_all_files)
    print(len(bins))
    print(len(plot_list_all_files))
    # bins = np.array(bins)

    # # fig = plt.figure()
    # # fig, ax = plt.subplots(1,1, figsize=(60,60))
    # # heatplot = ax.imshow(values_to_print, cmap='BuPu')
    # # ax.set_xticklabels(all_otus)
    # # ax.set_yticklabels(rbox_identifiers_all_drive)
    #
    #
    # fig = go.Figure(data=go.Heatmap(
    #     z=plot_list_all_files.transpose(),
    #     x=bins,
    #     y=rbox_protein_names,
    #     hoverongaps=False))
    # fig.update_layout(title_text="RBOX gene presence",
    #                   title_font_size=30)
    #
    # fig.update_layout(
    #     xaxis=dict(
    #         tickmode='linear',
    #         tick0=0.5,
    #         dtick=0.1
    #     )
    # )
    print(heatmap_annotations)
    # remove underscore from bacteria name and replace it with space
    for i, bin in enumerate(bins):
        bins[i] = bin.replace("_", " ")

    fig = ff.create_annotated_heatmap(z=plot_list_all_files, x=bins, y=protein_names,
                                      annotation_text=heatmap_annotations)

    # create annotation for heatmap 0 becomes - and 1 becomes X
    heatmap_annotation = []

    fig.show()


# call important functions to create heatmap sepcifying proteins to look for as well as directory for alignments
def compute_heatmap_for_directory(alignments_directory, protein_list_path, protein_names, alignment_file_name):
    bin_names = []
    bin_paths = []
    all_bins_rbox_proteins = []

    list_heatmap_all_bins = []
    annotation_heatmap_all_bins = []
    for (root, dirs, files) in os.walk(alignments_directory):
        if "name" in root and os.path.exists(root + alignment_file_name):
            bin_paths.append(root)
            bin = root.split("name")[1]
            bin_names.append(bin)

    for i, bin in enumerate(bin_names):
        try:
            all_bins_rbox_proteins.append(
                get_proteins_in_alignments(protein_list_path, bin_paths[i] + alignment_file_name, protein_names))

        except:
            print("File not found exception")

    for l in all_bins_rbox_proteins:
        list_heatmap_all_bins.append(create_list_for_plot(l, protein_names)[0])
        annotation_heatmap_all_bins.append(create_list_for_plot(l, protein_names)[1])

    create_heatmap(list_heatmap_all_bins, bin_names, annotation_heatmap_all_bins, protein_names)


# First create heatmap for all samples and all bins to get an overview of possible complete assemblies that have rbox proteins
for dir in glob.glob("/Users/timolucas/Documents/spirito/assembly/r*/extracted_reads/"):
    try:
        compute_heatmap_for_directory(dir,
                                      "/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                                      rbox_protein_names, "/r1t1_prokka.m8")
    except:
        print("Exception caught")

# then do heatmap again for combined directory with annotated bins that have high checkm completeness and low contamination
compute_heatmap_for_directory("/Users/timolucas/Documents/spirito/assembly/high_completeness/",
                              "/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                              rbox_protein_names, "/r1t1_prokka.m8")

# heatmaps for FAB proteins
# create heatmap for all of spiritos assemblies
# for dir in glob.glob("/Users/timolucas/Documents/spirito/assembly/r1t*/extracted_reads/"):
#   compute_heatmap_for_directory(dir,
#                                  "/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
#                                  fab_protein_names
#                                   , "/fab_alignments.m8")

# create heatmap for FAB proteins for bacteria with high completeness
compute_heatmap_for_directory("/Users/timolucas/Documents/spirito/assembly/high_completeness_fab/",
                              "/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
                              fab_protein_names, "/fab_alignments.m8")
