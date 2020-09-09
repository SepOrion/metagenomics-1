# This scripts takes a diamond .m8 alignment output in the following form:
# qseqid qlen sseqid slen qstart qend sstart send evalue bitscore score length pident positive nident
# it then filters the file for alignments with at least 50 percent positive matches and creates a new file

file = open('/Users/timolucas/Documents/spirito/new_annotation/rbox_proteins_id90_query-cover90.m8', 'r')
output_file = open(
    '/Users/timolucas/Documents/spirito/new_annotation/rbox_proteins_id90_query-cover90_positive_over50percent.m8', 'w')
for alignment in file:
    # get sequence length of alignment
    sequence_length = alignment.split('\t')[3]
    # if sequence length is not query length take longer sequence length as input for percent computation
    # get number of bases that have a positive alignment score for alignment
    positives = alignment.split('\t')[13]
    percent_positives = int(positives) / int(sequence_length)
    print(percent_positives)
    if (percent_positives > 0.5):
        print('asdf')
