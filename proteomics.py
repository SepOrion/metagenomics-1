import csv

from Bio import SeqIO

gene_ids_and_tax_R1 = []
gene_ids_and_tax_R2 = []
gene_ids_and_tax_R3 = []
gene_id_tax_dict = {}

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
                     "acetyl-CoA decarbonylase/synthase",
                     "acetyl-CoA decarbonylase",
                     "acetyl-CoA carboxylase",
                     "biotin carboxylase",
                     "S-malonyltransferase",
                     "3-oxoacyl-ACP synthase",
                     "3-oxoacyl-ACP reductase", "3-hydroxyacyl-ACP dehydratase"
    , "enoyl-ACP reductase",
                     "beta-ketoacyl-ACP synthase"]


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
        id = line_list[0]
        protein = line_list[2]
        alignment_length = line_list[11]
        alignment_num_positives = line_list[13]
        percent_positives = int(alignment_num_positives) / int(alignment_length)
        # only consider alignments with at least 50% positive aligned bases
        if percent_positives >= 0.5:
            functional_proteins_aligned.append((protein, id))

    # print(rbox_proteins_aligned)
    rbox_proteins_aligned_tax = []
    # get the names of the sequences that alignments were found for
    otu_proteins_dict = {}
    for record in SeqIO.parse(fasta_with_rbox_proteins, "fasta"):
        for protein in functional_proteins_aligned:
            if protein[0] == record.name:
                rbox_proteins_aligned_function.append(record.description)
                rbox_proteins_aligned_tax.append(protein[1])
    # deduplicate list
    # rbox_proteins_aligned_function = list(dict.fromkeys(rbox_proteins_aligned_function))
    rbox_proteins_aligned_trimmed = []
    # trim accessions

    # this function gets the .m8 alignment file and parses all

    rbox_proteins_with_alignments = []

    # check if rbox proteins are in clostridiales for alignment file and add them to new list
    for protein in considered_proteins:
        for i, rbox_protein_aligned in enumerate(rbox_proteins_aligned_function):
            if protein in rbox_protein_aligned:
                rbox_proteins_with_alignments.append((protein, gene_id_tax_dict[rbox_proteins_aligned_tax[i]]))

    # deduplicate list with rbox proteins
    # rbox_proteins_with_alignments = list(dict.fromkeys(rbox_proteins_with_alignments))
    rbox_proteins_with_alignments = list(dict.fromkeys(rbox_proteins_with_alignments))
    return (rbox_proteins_with_alignments)


def split_proteomics_fasta():
    with open(
            "/Users/timolucas/Google Drive/caprylate reactor paper/metagenomics and proteomics/Proteomics/AllSamplesCombinedFromSummer2016-Percolator/ProteinIDsSummer2016_csv.csv") as proteomics_csv:
        readCSV = csv.reader(proteomics_csv, delimiter=';')
        # parse csv file and split csv in three lists corresponding to reactor
        for row in readCSV:
            gene_id_tax_dict[row[0]] = row[16]
            if row[17] == "1":
                # add tuple of gene ID and taxonomic assignment to lists
                gene_ids_and_tax_R1.append((row[0], row[16]))
            if row[17] == "2":
                # add tuple of gene ID and taxonomic assignment to lists
                gene_ids_and_tax_R2.append((row[0], row[16]))
            if row[17] == "3":
                # add tuple of gene ID and taxonomic assignment to lists
                gene_ids_and_tax_R3.append((row[0], row[16]))
    # print(gene_ids_and_tax_R1)
    gene_id_only_R1 = ([lis[0] for lis in gene_ids_and_tax_R1])
    gene_id_only_R2 = ([lis[0] for lis in gene_ids_and_tax_R2])
    gene_id_only_R3 = ([lis[0] for lis in gene_ids_and_tax_R3])

    # next load fasta with protein sequences and split them according the reactors

    records_R1 = []
    records_R2 = []
    records_R3 = []

    for record in SeqIO.parse(
            "/Users/timolucas/Google Drive/caprylate reactor paper/metagenomics and proteomics/Proteomics/ProteinsWeIDedInProteomics.fasta",
            "fasta"):
        if record.id in gene_id_only_R1:
            records_R1.append(record)
        if record.id in gene_id_only_R2:
            records_R2.append(record)
        if record.id in gene_id_only_R3:
            records_R3.append(record)

    # write filtered fastas to .fasta files
    with open("/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R1.fasta", "w") as output_handle:
        SeqIO.write(records_R1, output_handle, "fasta")
    with open("/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R2.fasta", "w") as output_handle:
        SeqIO.write(records_R2, output_handle, "fasta")
    with open("/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R3.fasta", "w") as output_handle:
        SeqIO.write(records_R3, output_handle, "fasta")


split_proteomics_fasta()
# print(gene_id_tax_dict)
print("RBOX proteins in proteomics data: ")
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R1.m8",
                                 rbox_protein_names))
print(len(
    get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                               "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R1.m8",
                               rbox_protein_names)))
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R2.m8",
                                 rbox_protein_names))
print(len(
    get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                               "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R2.m8",
                               rbox_protein_names)))
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R3.m8",
                                 rbox_protein_names))
print(len(
    get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                               "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_rbox_R3.m8",
                               rbox_protein_names)))
print("FAB proteins in proteomics data: ")
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R1_FAB.m8",
                                 fab_protein_names))
print(len(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
                                     "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R1_FAB.m8",
                                     fab_protein_names)))
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R2_FAB.m8",
                                 fab_protein_names))
print(len(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/FAB/FAB_proteins.faa",
                                     "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R2_FAB.m8",
                                     fab_protein_names)))
print(get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                                 "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R3_FAB.m8",
                                 fab_protein_names))
print(len(
    get_proteins_in_alignments("/Users/timolucas/Documents/spirito/new_annotation/clear_rbox_proteins_filtered.fa",
                               "/Users/timolucas/Documents/spirito/proteomics/proteomics_seqs_R3_FAB.m8",
                               fab_protein_names)))
