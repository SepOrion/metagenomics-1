import glob

import pandas as pd
import plotly.express as px

for file in glob.glob("/Users/timolucas/Desktop/rbox_review_abundance_plots/*.daa.txt"):
    print(file)
    count_dict = {}
    read_count_file = open(file, 'r')
    for line in read_count_file.read().splitlines():
        taxonomy = line.split("\t")[1]
        count = line.split("\t")[2]
        count_dict[taxonomy] = count
    read_count_file.close()
    df = pd.DataFrame(count_dict.items(), columns=['taxonomy', 'count'])
    fig = px.pie(df, values='count', names='taxonomy', title='Abundance plot')
    fig.show()
