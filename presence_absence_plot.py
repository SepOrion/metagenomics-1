import os

kegg_dict = {}
kegg_enzyme_cur = []
otu_names_all = []
otu_names_cur = []
directory = os.fsencode(".")


def compute_dict_line(otu_names_cur, name, l):
    kegg_dict_temp = kegg_dict.copy()
    tuple_list = []
    for index, otu in enumerate(otu_names_cur):
        otu_val_tuple = (otu, int(l[index]))
        tuple_list.append(otu_val_tuple)
    kegg_dict_temp[name] = tuple_list
    return kegg_dict_temp


# filename = os.fsdecode(file)
filename = 'spirito_r1_t1-kegg.txt'


def compute_dict_file(filename):
    csv = open(filename)
    header = csv.readline()
    otu_names_cur = header.split('\t')
    otu_names_cur.pop(0)
    for line in csv.readlines():
        line_split = line.split('\t')
        name = line_split[0]
        # print(name)
        line_split.pop(0)
        dict = compute_dict_line(otu_names_cur, name, line_split)
        csv.close()
    return dict


dict_r1_t1 = compute_dict_file('spirito_r1_t1-kegg.txt')
dict_r1_t2 = compute_dict_file('spirito_r1_t2-kegg.txt')
dict_r1_t3 = compute_dict_file('spirito_r1_t3-kegg.txt')
dict_r2_t1 = compute_dict_file('spirito_r2_t1-kegg.txt')
dict_r2_t2 = compute_dict_file('spirito_r2_t2-kegg.txt')
dict_r2_t3 = compute_dict_file('spirito_r2_t3-kegg.txt')
dict_r3_t1 = compute_dict_file('spirito_r3_t1-kegg.txt')
dict_r3_t2 = compute_dict_file('spirito_r3_t2-kegg.txt')
dict_r3_t3 = compute_dict_file('spirito_r3_t3-kegg.txt')

dict_list = [dict_r1_t1, dict_r1_t2, dict_r1_t3, dict_r2_t1, dict_r2_t2, dict_r2_t3, dict_r3_t1, dict_r3_t2, dict_r3_t3]

print(dict_r1_t1)
