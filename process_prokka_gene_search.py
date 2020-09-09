import glob

proteins_from_prokka_dict = {}
output = open("/Users/timolucas/Documents/jeon/new_pathway_analysis/7_new_pathway/gene_search_sensitive.txt", 'w')

for prokka_filename in glob.glob("/Users/timolucas/Documents/jeon/new_pathway_analysis/7_new_pathway/*/*.faa"):
    if prokka_filename.endswith(".faa"):
        sample = prokka_filename.split("/")[7]
        sample = prokka_filename.split("/")[7]
        species = prokka_filename.split("/")[8][5:-4]
        alignment_filename = prokka_filename.replace(".faa", "_sensitive.m8")
        alignment_file = open(
            alignment_filename,
            'r')
        prokka_file = open(
            prokka_filename,
            'r')
        for line in prokka_file:
            if line.startswith(">"):
                annotated_protein_name = " ".join(line.split(" ")[1:])
                annotated_protein_id = line.split(" ")[0][1:]
                proteins_from_prokka_dict[annotated_protein_id] = annotated_protein_name

        for line in alignment_file:
            query = line.split("\t")[0]
            hit = line.split("\t")[2]
            print(sample + "" + species + " " + query + " " + hit + " " + proteins_from_prokka_dict[hit])
            output.write(sample + "\tt" + species + "\t" + query + "\t" + hit + "'\t" + proteins_from_prokka_dict[hit])
output.close()

# parse output again to count occurences of each gene
output = open("/Users/timolucas/Documents/jeon/new_pathway_analysis/7_new_pathway/gene_search_sensitive.txt", 'r')
