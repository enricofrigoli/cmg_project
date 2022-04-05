from os import listdir
import pandas as pd
import re
# from plotnine import ggplot, aes, geom_bar, position_dodge
import matplotlib.pyplot as plt

listToDf = list()

myPath = "./prokkaOutput"

for folder in listdir(myPath):
    # print("Folder: "+ folder)
    for filename in listdir("/".join([myPath,folder])):
        if filename.endswith("tsv"):
            cur = "/".join([myPath, folder, filename])
            curD = dict()
            with open(cur, "r") as c:
                curT = c.read()
            nCDS = len(re.findall("\tCDS\t", curT))
            nHyp = len(re.findall("hypothetical protein", curT))
            curD["DATA"] = filename[:-8]
            curD["CDS"] = nCDS
            curD["HYP"] = nHyp
            curD["PROT"] = nCDS - nHyp
            listToDf.append(curD)

df2plot = pd.DataFrame(listToDf)

print(df2plot)

plt.bar(df2plot.DATA, df2plot.PROT, color='orange')
plt.bar(df2plot.DATA, df2plot.HYP, bottom=df2plot.PROT, color='b')
plt.legend(["Known Proteins", "Hypothetical proteins"], loc="lower right")
plt.xlabel("Dataset")
plt.xlabel("MAGs")
plt.title("Number of CDS found per MAG")
plt.xticks(df2plot.DATA, rotation=90, fontsize=5)
plt.tight_layout()
plt.show()