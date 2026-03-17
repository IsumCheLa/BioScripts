import Bio
import os
from Bio import SeqIO
from tqdm.notebook import tqdm
from matplotlib import pyplot as plt

#ïåðåâîäèì êîëè÷åñòâî áàçîâûõ ïàð â êèëîáàéòû
def len2bytes(bp):
    return 2*bp/8/1024
    
#ñîçäàåì ìàññèâû, â êîòîðûå ïîëîæèì ïëîòíîñòè
rations_first_chrom = []
rations_second_chrom = []
#Ïåðåáèðàåì âñå îðãàíèçìû
for dir in tqdm( os.listdir( r"C:\Users\User\ÏÐÎÅÊÒÁÈÎ\ncbi_datasetREF2\ncbi_dataset\data")):
    path = os.path.join(r"C:\Users\User\ÏÐÎÅÊÒÁÈÎ\ncbi_datasetREF2\ncbi_dataset\data", dir, 'genomic.gbff')
    parsq = SeqIO.parse(path, 'genbank')
#íàñ èíòåðåñóþò òîëüêî îðãàíèìçû, ÷åé ãåíîì ñîñòîèò òîëüêî èç ïåðâè÷íîé è âòîðè÷íîé õðîìîñîì
#îñòàëüíîå ôèëüòðóåòñÿ
    fl = False
    for e, record in enumerate(parsq):
        if not(("chromosome I," in record.description) or ("chromosome II," in record.description) or ("chromosome 1," in record.description) or ("chromosome 2," in record.description)):
            fl = True


    if fl:
        continue
    if e != 1:
        continue 
    parsq = SeqIO.parse(path, 'genbank')
    for e, record in enumerate(parsq):
        #èùåì êîëè÷åñòâî ãåíîâ
        counter = 0
        for feature in record.features:
            if feature.type == "gene":
                counter += 1
            

        #äåëèì êîëè÷åñòâî ãåíîâ íà ðàçìåð õðîìîñîìû
        if ("chromosome II" in record.description) or ("chromosome 2" in record.description):
            rations_second_chrom.append(counter/len2bytes(len(record.seq)))
        else:
            rations_first_chrom.append(counter/len2bytes(len(record.seq)))
            

#Ðèñóåì ãðàôèê
plt.scatter(rations_first_chrom, rations_second_chrom, facecolors='none', edgecolors='grey')
plt.xlabel("chromosome I, gene number / chromosome size, Kb ")
plt.ylabel("chromosome II, gene number / chromosome size, Kb ")

x1, y1 = [3.65, 3.65], [4.2, 3.4]
x2, y2 = [3.65, 4.2], [3.4, 3.4]
plt.plot(x1, y1, x2, y2, linestyle = '--', color = "grey")
plt.savefig("plt2.png")
plt.show() #hm
