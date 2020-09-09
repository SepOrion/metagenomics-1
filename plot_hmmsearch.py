import glob

import Bio.SearchIO.HmmerIO as hmmer
import numpy as np
import pandas as pd
import plotly
import plotly.figure_factory as ff

files = sorted(glob.glob("/Users/timolucas/Documents/jeon/review/out_c6/*.txt"))

lps_domains = ["LptD", "LptF_LptG", "Prenyltransf", "LpxD", "Hexapep", "Hexapep_2", "Acetyltransf_11", "LpxB",
               "Sigma70_r4", "Sigma70_r4_2", "Sigma70_ECF", "Epimerase", "CTP_transf_3", "Glyco_trans_1_4",
               "Glycos_transf_1", "LpxC", "CBS", "SIS", "Hydrolase_3", "Hydrolase", "LptD_N", "ABC_tran", "EptA_B_N",
               "Sulfatase",
               "PagL", "Glycos_transf_N", "Glycos_transf_1", "Kdo", "Glycos_transf_N", "Glyco_transf_4",
               "Glyco_transf_9"]

pfam_domain_to_product_dict = {}

rbox_domains = ["3HCDH", "3HCDH_N", "Acyl-CoA_dh_1", "Acyl-CoA_dh_2", "Acyl-CoA_dh_M", "Acyl-CoA_dh_N", 'ECH_1',
                'ECH_2',
                'ETF_alpha', 'ETF', 'Thiolase_C', 'Thiolase_N', 'AcetylCoA_hyd_C', 'AcetylCoA_hydro']

pfam_domain_to_product_dict["LptD"] = "LPS-assembly protein LptD"
pfam_domain_to_product_dict["LptF_LptG"] = "Lipopolysaccharide export system permease protein LptF/LptG"
pfam_domain_to_product_dict[
    "Prenyltransf"] = "Ditrans,polycis-undecaprenyl-diphosphate synthase ((2E,6E)-farnesyl-diphosphate specific)"
pfam_domain_to_product_dict["LpxD"] = "UDP-3-O-(3-hydroxymyristoyl)glucosamine N-acyltransferase"
pfam_domain_to_product_dict[
    "Hexapep"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase / UDP-3-O-acylglucosamine N-acyltransferase 1"
pfam_domain_to_product_dict[
    "Hexapep_2"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase / UDP-3-O-acylglucosamine N-acyltransferase 1"
pfam_domain_to_product_dict[
    "Acetyltransf_11"] = "Acyl-[acyl-carrier-protein]--UDP-N-acetylglucosamine O-acyltransferase"
pfam_domain_to_product_dict["LpxB"] = "Lipid-A-disaccharide synthase"
pfam_domain_to_product_dict["Sigma70_r2"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_r4"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_r4_2"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Sigma70_ECF"] = "ECF RNA polymerase sigma-E factor"
pfam_domain_to_product_dict["Epimerase"] = "UDP-N-acetylglucosamine 4-epimerase"
pfam_domain_to_product_dict["Epimerase_2"] = "UDP-N-acetylglucosamine 4-epimerase"
pfam_domain_to_product_dict["CTP_transf_3"] = "3-deoxy-manno-octulosonate cytidylyltransferase"
pfam_domain_to_product_dict["Glyco_trans_1_4"] = "Glycogen(starch) synthase"
pfam_domain_to_product_dict["Glyco_transf_1"] = "Glycogen(starch) synthase"
pfam_domain_to_product_dict["LpxC"] = "UDP-3-O-acyl-N-acetylglucosamine deacetylase"
pfam_domain_to_product_dict["CBS"] = "Arabinose-5-phosphate isomerase"
pfam_domain_to_product_dict["SIS"] = "Arabinose-5-phosphate isomerase"
pfam_domain_to_product_dict["Hydrolase_3"] = "3-deoxy-manno-octulosonate-8-phosphatase"
pfam_domain_to_product_dict["Hydrolase"] = "3-deoxy-manno-octulosonate-8-phosphatase"
pfam_domain_to_product_dict["LptD_N"] = "Lipopolysaccharide export system protein LptA"
pfam_domain_to_product_dict["ABC_tran"] = "Hemin import ATP-binding protein HmuV"
pfam_domain_to_product_dict["EptA_B_N"] = "Phosphoethanolamine transferase CptA"
pfam_domain_to_product_dict["Sulfatase"] = "Phosphoethanolamine transferase CptA"
pfam_domain_to_product_dict["PagL"] = "Lipid A deacylase PagL"
pfam_domain_to_product_dict["Glycos_transf_N"] = "Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase"
pfam_domain_to_product_dict["Glycos_transf_1"] = "Lipid IV(A) 3-deoxy-D-manno-octulosonic acid transferase"
pfam_domain_to_product_dict["Kdo"] = "Lipopolysaccharide core heptose(I) kinase RfaP"
pfam_domain_to_product_dict["Glycos_transf1"] = "Lipopolysaccharide core biosynthesis protein RfaG"
pfam_domain_to_product_dict["Glyco_transf_4"] = "Lipopolysaccharide core biosynthesis protein RfaG"
pfam_domain_to_product_dict["Glyco_transf_9"] = "Lipopolysaccharide heptosyltransferase 1"

pfam_description = [pfam_domain_to_product_dict[x] for x in lps_domains]
# test = [pfam_description[pfam_description]]
print(pfam_description)
pfam_description = np.array(pfam_description)


# pfam_stacked
def get_query_results(hmmer_file):
    results = []
    file = open(hmmer_file, 'r')
    hmmer_parser = hmmer.hmmer3_domtab.Hmmer3TabParser(file)

    for query in hmmer_parser:
        print(query.id)
        print(len(query))
        print(query)
        results.append((query.id, len(query)))
    # results = list(set(results))
    print(results)

    return results


def get_otu_name_from_filename(file_name: str):
    otu_name = file_name.split("/")[-1]
    otu_name = otu_name.split(".")[0]
    otu_name = otu_name.split("name")[1]
    return otu_name


# create dictionary with OTU as keys and list with query hits as values
gene_presence_dict = {}


def get_results_and_plot(pfam_domain_list):
    sum_genes = 0
    heatmap_text = []
    for file in files:
        if pfam_domain_list == "pfam_domains":
            otu = get_otu_name_from_filename(file)
        else:
            otu = file.split("/")[-1]
            otu = otu.split(".faa")[0]
            # otu = otu.split("_")[1]
            # otu = otu.replace("-"," ")
        # fille dictionary with all hmmer search results in folder

        # get query results for all files and check if pfam id is present in query result - store result in list and then dict for plotting

        gene_presence_list = []
        heatmap_annotation_one_file = []
        results = get_query_results(file)
        if len(results) == 0:
            continue
        print("Results")
        print(results)
        query_result_names = list(list(zip(*results))[0])
        query_result_counts = list(list(zip(*results))[1])
        # print(query_result_names)
        # print(query_result_counts)

        for pdom in pfam_domain_list:
            for id, name in enumerate(query_result_names):
                if name == pdom:
                    gene_presence_list.append(query_result_counts[id])

            if pdom not in query_result_names:
                gene_presence_list.append(0)
        heatmap_text.append(gene_presence_list)

        print(heatmap_text)
        sum_genes += sum(gene_presence_list)
        # skip file if gene_presence list is empty = no query hits for for otu
        if sum(gene_presence_list) < 10:
            continue
        gene_presence_dict[otu] = gene_presence_list

    print(f"Total sum of gene counts: {sum_genes}")
    df = pd.DataFrame.from_dict(gene_presence_dict)
    df = df.to_numpy()

    fig = ff.create_annotated_heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()))
    # fig.add_trace(go.Heatmap(z=df, y=pfam_domain_list, x=list(gene_presence_dict.keys()),zmin=0,colorscale="blackbody"))
    fig.update_layout(title="RBOX pfam domains C6 producing reference genomes", xaxis=dict(side='bottom'), width=500,
                      height=1000)

    plotly.offline.plot(fig, filename='/Users/timolucas/Desktop/rbox_scarborough.html', auto_open=False)

    # fig.update_xaxes(side="top")
    fig.show()


get_results_and_plot(rbox_domains)

#
# ['aggrnyl', 'agsunset', 'algae', 'amp', 'armyrose', 'balance',
#  'blackbody', 'bluered', 'blues', 'blugrn', 'bluyl', 'brbg',
#  'brwnyl', 'bugn', 'bupu', 'burg', 'burgyl', 'cividis', 'curl',
#  'darkmint', 'deep', 'delta', 'dense', 'earth', 'edge', 'electric',
#  'emrld', 'fall', 'geyser', 'gnbu', 'gray', 'greens', 'greys',
#  'haline', 'hot', 'hsv', 'ice', 'icefire', 'inferno', 'jet',
#  'magenta', 'magma', 'matter', 'mint', 'mrybm', 'mygbm', 'oranges',
#  'orrd', 'oryel', 'peach', 'phase', 'picnic', 'pinkyl', 'piyg',
#  'plasma', 'plotly3', 'portland', 'prgn', 'pubu', 'pubugn', 'puor',
#  'purd', 'purp', 'purples', 'purpor', 'rainbow', 'rdbu', 'rdgy',
#  'rdpu', 'rdylbu', 'rdylgn', 'redor', 'reds', 'solar', 'spectral',
#  'speed', 'sunset', 'sunsetdark', 'teal', 'tealgrn', 'tealrose',
#  'tempo', 'temps', 'thermal', 'tropic', 'turbid', 'twilight',
#  'viridis', 'ylgn', 'ylgnbu', 'ylorbr', 'ylorrd'].
