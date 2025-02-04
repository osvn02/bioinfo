from Bio import SeqIO
from Bio import Entrez
from Bio import motifs

with open("data/pfm.txt") as handle:
    retour = motifs.parse(handle, "jaspar")
# print(len(retour))

sequence_amont = SeqIO.read("./data/seq_amont.fa", "fasta")

"""
for matrice in retour:
    matrice_base = matrice.counts
    matrice_corrigee = matrice_base.normalize(0.1)
    matrice_pssm = matrice_corrigee.log_odds()
    print(matrice_pssm.calculate(sequence_amont.seq))
"""


matrice = retour[0].counts
print(type(matrice))
print(matrice)

matrice_corrigee = matrice.normalize(0.1)
print(type(matrice_corrigee))
print(matrice_corrigee)

matrice_pssm = matrice_corrigee.log_odds()
print(type(matrice_pssm))
print(matrice_pssm)

# print(matrice_pssm.calculate(sequence_amont.seq))
