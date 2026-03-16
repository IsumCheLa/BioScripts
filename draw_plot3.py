import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

df = pd.read_csv("chrmb.proteinortho.tsv", sep='\t')

counts, bins = np.histogram(df["# Species"].values, bins = 50)
plt.stairs(np.log(counts), bins, fill=True, color='grey')
plt.ylabel("log(count)")
plt.xlabel("number of species in row")
plt.savefig("hist.png")