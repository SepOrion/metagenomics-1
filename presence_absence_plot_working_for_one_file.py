import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go

kegg_dict = {}
kegg_enzyme_cur = []
otu_names_all = []
otu_names_cur = []
directory = os.fsencode(".")


def compute_dict_line(kegg_dict, otu_names_cur, name, l):
    tuple_list = []
    for index, otu in enumerate(otu_names_cur):
        otu_val_tuple = (otu, int(l[index]))
        tuple_list.append(otu_val_tuple)
    kegg_dict[name] = tuple_list


# filename = os.fsdecode(file)


def compute_dict_file(filename):
    csv = open(filename)
    header = csv.readline()
    otu_names_cur = header.split('\t')
    otu_names_all.append(otu_names_cur)
    otu_names_cur.pop(0)
    for line in csv.readlines():
        line_split = line.split('\t')
        name = line_split[0]
        # print(name)
        line_split.pop(0)
        compute_dict_line(kegg_dict, otu_names_cur, name, line_split)
    return kegg_dict.copy()


def search_dict_by_kegg_ident(kegg_identifier, dictionary):
    for key in dictionary.keys():
        if kegg_identifier in key:
            return dictionary[key]


# returns all bacteria with a given  identifier in dictionary
def search_bacteria_with_kegg_identifier(sum_dict, kegg_identifier):
    bacteria_with_identifier = []
    print("The following bacteria had reads mapped to the gene " + kegg_identifier + ":")
    for otu in search_dict_by_kegg_ident(kegg_identifier, sum_dict):
        if otu[1] > 0:
            print(otu[0])
            bacteria_with_identifier.append(otu[0])
    return bacteria_with_identifier


# prints all bacteria with
def search_bacteria_with_rbox(sum_dict, identifiers):
    bacteria_with_rbox_idents = []
    for ident in identifiers:
        bacteria_with_ident = search_bacteria_with_kegg_identifier(sum_dict, ident)
        if len(bacteria_with_ident) != 0:
            for bacterium in bacteria_with_ident:
                bacteria_with_rbox_idents.append(bacterium)

    for bacterium in bacteria_with_rbox_idents:
        bacteria_count_with_rbox = bacteria_with_rbox_idents.count(bacterium)
        if bacteria_count_with_rbox == len(identifiers):
            print("The bacterium " + bacterium + "has reads mapped to all the important RBOX genes.")
        if bacteria_count_with_rbox < len(identifiers) and (
                len(identifiers) - bacteria_count_with_rbox <= len(identifiers) / 2):
            print("The bacterium " + bacterium + "has reads mapped to at least half of the RBOX genes.")


def search_kegg_enzymes_for_bacterium(bacterium_name, sum_dict):
    temp_list = []
    kegg_enzymes_list = []
    for key in sum_dict.keys():
        temp_list.append((key, search_dict_by_kegg_ident(key, sum_dict)))
    for enzyme in temp_list:
        for tuple in enzyme[1]:
            if bacterium_name == tuple[0] and tuple[1] > 0:
                kegg_enzymes_list.append(enzyme[0])
    return kegg_enzymes_list


###makes it easier to search dicts using only the kegg identifier

print('Parsing kegg tables ...')
dict_r1_t1 = compute_dict_file('spirito_r1_t1-kegg.txt')
dict_r1_t2 = compute_dict_file('spirito_r1_t2-kegg.txt')
dict_r1_t3 = compute_dict_file('spirito_r1_t3-kegg.txt')
dict_r2_t1 = compute_dict_file('spirito_r2_t1-kegg.txt')
dict_r2_t2 = compute_dict_file('spirito_r2_t2-kegg.txt')
dict_r2_t3 = compute_dict_file('spirito_r2_t3-kegg.txt')
dict_r3_t1 = compute_dict_file('spirito_r3_t1-kegg.txt')
dict_r3_t2 = compute_dict_file('spirito_r3_t2-kegg.txt')
dict_r3_t3 = compute_dict_file('spirito_r3_t3-kegg.txt')
dict_proteomics = compute_dict_file('/Users/timolucas/Documents/spirito/kegg/proteomics_reads_mpi60_kegg_taxonname.txt')

# Dict list is the list that is in the end used for plotting
# to change the data used for plotting load the correct text file using the compute_dict_file function and add it to dict list

dict_list = [dict_r1_t1, dict_r1_t2, dict_r1_t3, dict_r2_t1, dict_r2_t2, dict_r2_t3, dict_r3_t1, dict_r3_t2, dict_r3_t3]
# dict_list = [dict_r1_t1, dict_r1_t2, dict_r1_t3, dict_r2_t1, dict_r2_t2, dict_r2_t3, dict_r3_t1, dict_r3_t2, dict_r3_t3]
# dict_list = [dict_r1_t1, dict_r1_t2, dict_r1_t3]
# dict_list = [dict_r2_t1, dict_r2_t2, dict_r2_t3]
# dict_list = [dict_r3_t1,dict_r3_t2,dict_r3_t3]
# dict_list = [dict_proteomics]

all_otus = []

# Get all otus and remove duplicates
all_otus = []
all_keys_temp = []
all_keys = []
for d in dict_list:
    all_keys_temp.append(list(d.keys()))
for key_list in all_keys_temp:
    for key in key_list:
        if key not in all_keys:
            all_keys.append(key)

print('')
for l in otu_names_all:
    for item in l:
        all_otus.append(item)
all_otus = list(dict.fromkeys(all_otus))
print(all_otus)

# Overwrite all_otus to only plot most abundant bacteria in presence absence plot ( at least 10,000 reads in 5 of the 9 samples)
all_otus = ["Bacteroides reticulotermitis", "Bacteroidaceae bacterium HV4-6-C5C", "Prevotella oryzae",
            "Prevotella paludivivens", "Desulfovibrio legallii", "Pseudoclavibacter soli", "Lactobacillus fuchuensis",
            "Lactobacillus sakei", "Clostridium amylolyticum", "Clostridium kluyveri", "Clostridium sp. JN500901",
            "Clostridium sp. KNHs216", "Clostridium sp. W14A", "Oscillibacter ruminantium",
            "Caproiciproducens galactitolivorans", "Caproiciproducens sp. NJN-50",
            "Pseudoflavonifractor sp. Marseille-P3106", "Ruminococcaceae bacterium HV4-5-B5C",
            "Clostridiales bacterium", "Methanobacterium paludis"]

sum_dict = {}

# init sum dict that contains summed up data of all tables

# initialize sum_dict keys
for k in all_keys:
    sum_dict[k] = list()
# init with all keys
for index, otu in enumerate(all_otus):
    for k in sum_dict.keys():
        sum_dict[k].append((otu, 0))
# sum up all values for the same otus and  keys
print('Merging kegg tables ...')
for d in dict_list:
    for key in sum_dict.keys():
        if key in d.keys():
            for i, otu in enumerate(sum_dict[key]):
                for otu_d in d[key]:
                    if otu[0] == otu_d[0]:
                        # print('a')
                        sum_dict[key][i] = (otu_d[0], (otu[1] + otu_d[1]))

# print(search_dict_by_kegg_ident('K00012', sum_dict))
# print(search_dict_by_kegg_ident('K00036',sum_dict))
# print(search_dict_by_kegg_ident('K00111',sum_dict))

# bacteria_k00012 = search_bacteria_with_kegg_identifier(sum_dict,"K00012")

# ----- PUT RBOX IDENTIFIERS IN THE rbox_identifiers list then call search_bacteria_with_rbox function with the list----

rbox_identifiers_jeff = ["K00001", "K04072", "K00626", "K07516", "K00074", "K17865", "K17865", "K00248", "K01034",
                         "K01073", "K00625", "K00925"]
rbox_identifiers_all_drive = ["K00001 alcohol dehydrogenase [EC:1.1.1.1]",
                              "K04072 acetaldehyde dehydrogenase / alcohol dehydrogenase [EC:1.2.1.10 1.1.1.1]",
                              "K00626 acetyl-CoA C-acetyltransferase [EC:2.3.1.9]",
                              "K07516 3-hydroxyacyl-CoA dehydrogenase [EC:1.1.1.35]",
                              "K00074 3-hydroxybutyryl-CoA dehydrogenase [EC:1.1.1.157]",
                              "K01782 3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase [EC:1.1.1.35 4.2.1.17 5.1.2.3]",
                              "K01782 3-hydroxyacyl-CoA dehydrogenase / enoyl-CoA hydratase / 3-hydroxybutyryl-CoA epimerase [EC:1.1.1.35 4.2.1.17 5.1.2.3]",
                              "K01692 enoyl-CoA hydratase [EC:4.2.1.17]",
                              "K17865 3-hydroxybutyryl-CoA dehydratase [EC:4.2.1.55]",
                              "K00249 acyl-CoA dehydrogenase [EC:1.3.8.7]",
                              "K00248 butyryl-CoA dehydrogenase [EC:1.3.8.1]",
                              "K01034 acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]",
                              "K01034 acetate CoA/acetoacetate CoA-transferase alpha subunit [EC:2.8.3.8 2.8.3.9]",
                              "K01073 acyl-CoA hydrolase [EC:3.1.2.20]",
                              "K01071 medium-chain acyl-[acyl-carrier-protein] hydrolase [EC:3.1.2.21]",
                              "K00625 phosphate acetyltransferase [EC:2.3.1.8]", "K00925 acetate kinase [EC:2.7.2.1]",
                              "K03614 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit D [EC:7.2.1.2]",
                              "K03615 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit C [EC:7.2.1.2]",
                              "K03616 Na+-translocating ferredoxin:NAD+ oxidoreductase subunit B [EC:7.2.1.2]",
                              "K00532 ferredoxin hydrogenase [EC:1.12.7.2]", "K17998",
                              "K03522 electron transfer flavoprotein alpha subunit",
                              "K03521 electron transfer flavoprotein beta subunit", "K14086 ech hydrogenase subunit A",
                              "K14087 ech hydrogenase subunit B"]

# rbox_identifiers_all_drive = ["K00001", "K04072", "K00626", "K07516", "K00074", "K01782", "K01782", "K01782",
#                                     "K01692",
#                                     "K17865", "K00249", "K00248", "K00209", "K00209", "K01034", "K01034", "K01073",
#                                     "K01071",
#                                     "K00625", "K00925", "K03614", "K03615", "K03616", "K00532", "K17998"]


# search_bacteria_with_rbox(sum_dict, rbox_identifiers)

# you can call the search_kegg_enzymes_for_bacterium function ( needs bacterium name and sum_dict as argument)
# it gives you back a list with
# kegg_enzymes_valericigenes = search_kegg_enzymes_for_bacterium('Oscillibacter valericigenes', sum_dict)
# kegg_enzymes_clostridium_kluyveri = search_kegg_enzymes_for_bacterium('Clostridium kluyveri', sum_dict)

# print("Kluyveri")  # 962 in total
# print(kegg_enzymes_clostridium_kluyveri)
# print(len(kegg_enzymes_clostridium_kluyveri))
# print("Valericigenes")  # 1065 in total
# print(kegg_enzymes_valericigenes)
# print(len(kegg_enzymes_valericigenes))


# extract alignment count from table for rbox identifiers for plotting
values_rbox_identifiers = []
for ident in rbox_identifiers_all_drive:
    print(ident)
    temp = search_dict_by_kegg_ident(ident, sum_dict)
    print(temp)
    # with this if condition it doesn't add kegg enzymes to the plot that are not in the table
    if temp is not None:
        values_rbox_identifiers.append(temp)

values_to_print = []
print_dict = {}

# heatmap with complete alignment count
# for val in values_rbox_identifiers:
#     temp = []
#     for i in val:
#         temp.append(i[1])
#     values_to_print.append(temp)

# heatmap with 0 and 1
for val in values_rbox_identifiers:
    temp = []
    for i in val:
        if i[1] > 0:
            temp.append(1)
        else:
            temp.append(0)
    values_to_print.append(temp)

print('')
# 22 enzymes and 98 otus = 2156 points to plot

values_to_print = np.array(values_to_print)
values_to_print = pd.DataFrame(values_to_print)

# fig = plt.figure()
# fig, ax = plt.subplots(1,1, figsize=(60,60))
# heatplot = ax.imshow(values_to_print, cmap='BuPu')
# ax.set_xticklabels(all_otus)
# ax.set_yticklabels(rbox_identifiers_all_drive)


fig = go.Figure(data=go.Heatmap(
    z=values_to_print,
    x=all_otus,
    y=rbox_identifiers_all_drive,
    hoverongaps=False))
fig.show()

# clear values_to_print for second plot (not binary)
values_to_print = []

for val in values_rbox_identifiers:
    temp = []
    for i in val:
        temp.append(i[1])
    values_to_print.append(temp)

fig = go.Figure(data=go.Heatmap(
    z=values_to_print,
    x=all_otus,
    y=rbox_identifiers_all_drive,
    hoverongaps=False))
fig.show()
