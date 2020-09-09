# This script takes a a FASTA file with protein reads from ncbi as input and creates a filtered output by providing a list of search terms

import os

from Bio import SeqIO

protein_list_THIO = ["acetyl-CoA C-acetyltransferase", "acetyl-CoA acetyltransferase", "acetyl-CoA acyltransferase",
                     "acetyl-CoA C-acyltransferase"]
protein_list_THIO1 = ["acetyl-CoA acyltransferase", "acetyl-CoA C-acyltransferase"]
protein_list_THIO2 = ["acetyl-CoA C-acetyltransferase", "acetyl-CoA acetyltransferase"]

protein_list_RBOX = ['alcohol dehydrogenase',
                     "Acetaldehyde dehydrogenase",
                     "3-Hydroxyacyl-CoA dehydrogenase",
                     "3-Hydroxybutyryl-CoA",
                     "acetyl-CoA C-acyltransferase",
                     "Beta-ketothiolase",
                     "3-Hydroxyacyl-CoA dehydrogenase",
                     "3-hydroxybutyryl-CoA epimerase",
                     "enoyl-CoA dehydratase",
                     "3-Hydroxybutyryl-CoA dehydratase"
                     "Acyl-CoA dehydrogenase",
                     "short-chain acyl-CoA dehydrogenase",
                     "Acyl-CoA dehydrogenase",
                     "short-chain acyl-CoA dehydrogenase",
                     "Acyl-CoA dehydrogenase",
                     "short-chain acyl-CoA dehydrogenase",
                     "Acyl-CoA dehydrogenase"
                     "short-chain acyl-CoA dehydrogenase",
                     "short-chain acyl-CoA dehydrogenase",
                     "Acyl-CoA dehydrogenase",
                     "acyl-CoA dehydrogenase",
                     "trans-2-enoyl-CoA reductase",
                     "enoyl-[acyl-carrier-protein] reductase (NADH)",
                     "enoyl-[acyl-carrier-protein] reductase",
                     "Acetate-CoA transferase",
                     "butyrate---acetoacetate CoA-transferase",
                     "Acyl-CoA hydrolase",
                     "medium-chain acyl-[acyl-carrier-protein] hydrolase",
                     "phosphate acetyltransferase",
                     "acetate kinase",
                     "Rnf - ferredoxin–NAD+ oxidoreductase–Na+ translocating",
                     "H2ase ferredoxin hydrogenase",
                     "Rnf - ferredoxin–NAD+ oxidoreductase–Na+ translocating",
                     "acetyl-CoA synthase",
                     "confurcating Nfn",
                     "alcohol dehydrogenase",
                     "aldehyde reductase",
                     "ADH",
                     "alcohol dehydrogenase (NAD)",
                     "aliphatic alcohol dehydrogenase",
                     "ethanol dehydrogenase",
                     "NAD-dependent alcohol dehydrogenase",
                     "NAD-specific aromatic alcohol dehydrogenase",
                     "NADH-alcohol dehydrogenase",
                     "NADH-aldehyde dehydrogenase",
                     "primary alcohol dehydrogenase",
                     "yeast alcohol dehydrogenase",
                     "acetaldehyde dehydrogenase (acetylating)",
                     "aldehyde dehydrogenase (acylating)",
                     "acylating acetaldehyde dehyrogenase",
                     "DmpF", "BphJ", "3-hydroxyacyl-CoA dehydrogenase",
                     "beta-hydroxyacyl dehydrogenase",
                     "beta-keto-reductase", "3-keto reductase", "3-hydroxyacyl coenzyme A dehydrogenase",
                     "beta-hydroxyacyl-coenzyme A synthetase", "beta-hydroxyacylcoenzyme A dehydrogenase",
                     "beta-hydroxybutyrylcoenzyme A dehydrogenase", "3-hydroxyacetyl-coenzyme A dehydrogenase",
                     "L-3-hydroxyacyl coenzyme A dehydrogenase", "L-3-hydroxyacyl CoA dehydrogenase",
                     "beta-hydroxyacyl CoA dehydrogenase",
                     "3beta-hydroxyacyl coenzyme A dehydrogenase", "3-hydroxybutyryl-CoA dehydrogenase",
                     "beta-ketoacyl-CoA reductase",
                     "beta-hydroxy acid dehydrogenase", "3-L-hydroxyacyl-CoA dehydrogenase",
                     "3-hydroxyisobutyryl-CoA dehydrogenase", "1-specific DPN-linked beta-hydroxybutyric dehydrogenase",
                     "3-hydroxybutyryl-CoA dehydrogenase", "beta-hydroxybutyryl coenzyme A dehydrogenase",
                     "L(+)-3-hydroxybutyryl-CoA dehydrogenase",
                     "BHBD;", "L-3-hydroxybutyryl coenzyme A (nicotinamide adenine dinucleotide phosphate)",
                     "L-(+)-3-hydroxybutyryl-CoA dehydrogenase", "beta-hydroxybutyryl-CoA dehydrogenase",
                     "acetyl-CoA C-acetyltransferase", "acetoacetyl-CoA thiolase",
                     "beta-acetoacetyl coenzyme A thiolase;",
                     "2-methylacetoacetyl-CoA thiolase [misleading]", "3-oxothiolase", "acetyl coenzyme A thiolase",
                     "acetyl-CoA acetyltransferase",
                     "acetyl-CoA:N-acetyltransferase", "thiolase II", "type II thiolase",
                     "acetyl-CoA C-acyltransferase",
                     "beta-ketothiolase", "3-ketoacyl-CoA thiolase", "KAT", "beta-ketoacyl coenzyme A thiolase",
                     "beta-ketoacyl-CoA thiolase",
                     "beta-ketoadipyl coenzyme A thiolase", "beta-ketoadipyl-CoA thiolase", "3-ketoacyl CoA thiolase",
                     "3-ketoacyl coenzyme A thiolase", "3-ketoacyl thiolase",
                     "3-ketothiolase", "3-oxoacyl-CoA thiolase", "3-oxoacyl-coenzyme A thiolase",
                     "6-oxoacyl-CoA thiolase",
                     "acetoacetyl-CoA beta-ketothiolase", "acetyl-CoA acyltransferase",
                     "ketoacyl-CoA acyltransferase", "ketoacyl-coenzyme A thiolase",
                     "long-chain 3-oxoacyl-CoA thiolase",
                     "oxoacyl-coenzyme A thiolase", "pro-3-ketoacyl-CoA thiolase", "thiolase I",
                     "type I thiolase", "2-methylacetoacetyl-CoA thiolase [misleading]",
                     "3-hydroxybutyryl-CoA epimerase",
                     "3-hydroxybutyryl coenzyme A epimerase",
                     "3-hydroxyacyl-CoA epimerase", "enoyl-CoA hydratase", "enoyl hydrase",
                     "unsaturated acyl-CoA hydratase",
                     "beta-hydroxyacyl-CoA dehydrase", "beta-hydroxyacid dehydrase",
                     "acyl coenzyme A hydrase", "crotonase", "crotonyl hydrase", "2-octenoyl coenzyme A hydrase",
                     "enoyl coenzyme A hydratase", "2-enoyl-CoA hydratase",
                     "short-chain enoyl-CoA hydratase", "ECH", "trans-2-enoyl-CoA hydratase",
                     "enoyl coenzyme A hydrase (D)",
                     "enoyl coenzyme A hydrase (L)",
                     "short chain enoyl coenzyme A hydratase", "D-3-hydroxyacyl-CoA dehydratase", "enol-CoA hydratase",
                     "3-hydroxybutyryl-CoA dehydratase",
                     "D-3-hydroxybutyryl coenzyme A dehydratase", "D-3-hydroxybutyryl-CoA dehydratase",
                     "enoyl coenzyme A hydrase(D)",
                     "(3R)-3-hydroxybutanoyl-CoA hydro-lyase", "3-hydroxybutyryl-CoA dehydratase",
                     "D-3-hydroxybutyryl coenzyme A dehydratase",
                     "D-3-hydroxybutyryl-CoA dehydratase",
                     "enoyl coenzyme A hydrase (D)",
                     "(3R)-3-hydroxybutanoyl-CoA hydro-lyase", "short-chain acyl-CoA dehydrogenase",
                     "butyryl-CoA dehydrogenase",
                     "butanoyl-CoA dehydrogenase",
                     "butyryl dehydrogenase",
                     "unsaturated acyl-CoA reductase",
                     "ethylene reductase",
                     "enoyl-coenzyme A reductase",
                     "unsaturated acyl coenzyme A reductase",
                     "butyryl coenzyme A dehydrogenase",
                     "short-chain acyl CoA dehydrogenase",
                     "short-chain acyl-coenzyme A dehydrogenase",
                     "3-hydroxyacyl CoA reductase",
                     "butanoyl-CoA:(acceptor) 2,3-oxidoreductase",
                     "ACADS (gene name)", "enoyl-[acyl-carrier-protein] reductase (NADH)",
                     "enoyl-[acyl carrier protein] reductase",
                     "enoyl-ACP reductase",
                     "NADH-enoyl acyl carrier protein reductase",
                     "NADH-specific enoyl-ACP reductase",
                     "acyl-[acyl-carrier-protein]:NAD+ oxidoreductase",
                     "fabI (gene name",
                     "acetate CoA-transferase",
                     "acetate coenzyme A-transferase",
                     "butyryl CoA:acetate CoA transferase",
                     "butyryl coenzyme A transferase",
                     "butyrate---acetoacetate CoA-transferase",
                     "butyryl coenzyme A-acetoacetate coenzyme A-transferase",
                     "butyryl-CoA-acetoacetate CoA-transferase",
                     "acetyl-CoA hydrolase",
                     "acetyl-CoA deacylase",
                     "acetyl-CoA acylase",
                     "acetyl coenzyme A hydrolase",
                     "acetyl coenzyme A deacylase",
                     "acetyl coenzyme A acylase",
                     "acetyl-CoA thiol esterase",
                     "dodecanoyl-[acyl-carrier-protein] hydrolase",
                     "lauryl-acyl-carrier-protein hydrolase",
                     "dodecanoyl-acyl-carrier-protein hydrolase",
                     "dodecyl-acyl-carrier protein hydrolase",
                     "dodecanoyl-[acyl-carrier protein] hydrolase",
                     "dodecanoyl-[acyl-carrier-protein] hydrolase",
                     "phosphate acetyltransferase",
                     "phosphotransacetylase",
                     "phosphoacylase",
                     "acetate kinase",
                     "acetokinase",
                     "acetic kinase",
                     "acetate kinase (phosphorylating)",
                     "acetate---CoA ligase",
                     "acetyl-CoA synthetase",
                     "acetyl activating enzyme",
                     "acetate thiokinase",
                     "acyl-activating enzyme",
                     "acetyl coenzyme A synthetase",
                     "acetic thiokinase",
                     "acetyl CoA ligase",
                     "acetyl CoA synthase",
                     "acetyl-coenzyme A synthase",
                     "short chain fatty acyl-CoA synthetase",
                     "short-chain acyl-coenzyme A synthetase",
                     "ferredoxin---NAD+ oxidoreductase (Na+-transporting)",
                     "Rnf complex (ambiguous)",
                     "Na+-translocating ferredoxin:NAD+ oxidoreductase",
                     "ferredoxin hydrogenase",
                     "H2 oxidizing hydrogenase",
                     "H2 producing hydrogenase [ambiguous]",
                     "NAD(P)+ transhydrogenase (ferredoxin)",
                     "NADH-dependent reduced ferredoxin:NADP+ oxidoreductase",
                     "Nfn",
                     "NfnAB",
                     "nfnAB (gene names)", "electron transfer flavoprotein", "etfA", "etfB",
                     "electron transfer flavoprotein subunit"]

protein_list_FAB = ["alcohol dehydrogenase", "aldehyde reductase", "ethanol dehydrogenase",
                    "anaerobic carbon-monoxide dehydrogenase", "carbon-monoxide dehydrogenase"
                                                               "Acetyl-CoA decarbonylase/synthase",
                    "Acetyl-CoA decarbonylase",
                    "acetyl-CoA carboxylase",
                    "biotin carboxylase",
                    "[acyl-carrier-protein] S-malonyl transferase", "S-malonyl transferase", "malonyltransferase", ""
                                                                                                                   "3-oxoacyl-[acyl-carrier-protein] synthase II",
                    "3-oxoacyl-[acyl-carrier-protein] synthase Ⅲ"
    , "3-oxoacyl-[acyl-carrier-protein] reductase", "3-hydroxyacyl-[acyl-carrier-protein] dehydratase"
    , "enoyl-[acyl-carrier-protein] reductase", "malonyl transferase", "S-malonyl transferase", "malonyl transacylase",
                    "biotin-carboxyl-carrier-protein:carbon-dioxide ligase",
                    "beta-ketoacyl-[acyl-carrier-protein] synthase II", "KASII", "KAS II",
                    "D-3-hydroxyoctanoyl-[acyl carrier protein] dehydratase", "D-3 dehydratase", "D-3-ACP dehydratase",
                    "3-hydroxyacyl-ACP dehydratase", "3-oxoacyl-ACP synthase"]


# iterate over all FASTA files in directory


def sequence_cleaner(fasta_file, min_length=0, por_n=100):
    # Create our hash table to add the sequences
    sequences = {}

    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        if (len(sequence) >= min_length and
                (float(sequence.count("N")) / float(len(sequence))) * 100 <= por_n):
            # If the sequence passed in the test "is it clean?" and it isn't in the
            # hash table, the sequence and its id are going to be in the hash
            if sequence not in sequences:
                sequences[sequence] = seq_record.id + "_" + seq_record.description
            # If it is already in the hash table, we're just gonna concatenate the ID
            # of the current sequence to another one that is already in the hash table
            else:
                sequences[sequence] += "_" + seq_record.id

    # Write the clean sequences

    # Create a file in the same directory where you ran this script
    with open(fasta_file + "_clear", "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + sequences[sequence] + "\n" + sequence + "\n")

    print("CLEAN!!!\nPlease check clear_" + fasta_file)


# given a list of protein names this function filters all fasta files in the specified directory for the proteins in the list and creates filtered fasta files
def filter_fastas_for_proteins(filter_list, filtered_fasta_suffix, directory):
    for filename in os.listdir(directory):
        if (filename.endswith(".faa") or filename.endswith(".fa") or filename.endswith(
                ".fasta")) and ("_filtered_" not in filename and "RBOX" not in filename):
            print(directory, os.path.join(filename))
            fasta_name = filename
            filtered = []

            # For loop over fasta file - adds all record of fasta to filtered list filtering by protein_list
            for record in SeqIO.parse(
                    directory + fasta_name, "fasta"):
                for protein in filter_list:
                    if protein.lower() in record.description.lower():
                        filtered.append(record)

            out_name = directory + fasta_name + '_filtered_' + filtered_fasta_suffix + ".faa"
            SeqIO.write(filtered, out_name, "fasta")
    return out_name


def filter_prokka_annotations(directory):
    dir_files = os.listdir(directory)
    for bin_folder in dir_files:
        if "read" in bin_folder and ".fasta" not in bin_folder:
            filtered_file_path = directory + bin_folder + "/"
            filtered_name = filter_fastas_for_proteins(protein_list_RBOX, "_RBOX", filtered_file_path)
            sequence_cleaner(filtered_name)
            filtered_name = filter_fastas_for_proteins(protein_list_THIO, "_THIOLASES", filtered_file_path)
            sequence_cleaner(filtered_name)


# filter .faa annotations from prokka for rbox genes for a given directory containing directories with annotatd bins created with MEGAN read exfractor

# filter_fastas_for_proteins(protein_list_RBOX,"rbox",'/Users/timolucas/Documents/spirito/new_annotation/')
# filter_fastas_for_proteins(protein_list_FAB,"fab",'/Users/timolucas/Documents/spirito/new_annotation/')
# filter_fastas_for_proteins(protein_list_THIO1,"thio1",'/Users/timolucas/Documents/spirito/new_annotation/')

# filter_fastas_for_proteins(protein_list_THIO1,"thio1",'/Users/timolucas/Documents/spirito/new_annotation/')


# filter prokka annotations folders for RBOX proteins
directories = ["/Users/timolucas/Documents/jeon/2n3_flye_medaka_extracted_reads/",
               "/Users/timolucas/Documents/jeon/4n5n6_flye_medaka_extracted_reads/",
               "/Users/timolucas/Documents/jeon/7n8_flye_medaka_extracted_reads/"]

# filter FAB proteins from NCBI


# for directory in directories:
#      filter_prokka_annotations(directory)


filter_fastas_for_proteins(protein_list_FAB, "FAB", "/Users/timolucas/Documents/spirito/new_annotation/FAB/")
