import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("/Users/timolucas/Documents/jeon/2n3_abundances.csv")
print(data.head())

data.plot(kind='bar', stacked=True)
plt.show()
