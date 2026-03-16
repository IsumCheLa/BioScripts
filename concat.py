from Bio import AlignIO
from Bio import SeqIO

#загружаем выравнивания
aligns = []
for f in os.listdir("./ort_align/"):
    aligns.append(AlignIO.read(os.path.join("./ort_align/", f), "fasta"))

#приводим к тому чтобы имена шли в одном порядке
for i in range(len(aligns)):
    aligns[i].sort()

#конкатинируем
s = aligns[0]
for i in range(1, len(aligns)):
    s += aligns[i]
#сохраняем
count = SeqIO.write(s, "example.afa", "fasta")    

#raxml не любит эти символы удалим их.
with open("example.afa", "r") as f:
    txt = f.read()

change = [ ":", ",", ")", "(", ";", "]", "[", "'"]

for sm in change:
    txt = txt.replace(sm, '_')

with open("example3.afa", "w") as f:
    f.write(txt)