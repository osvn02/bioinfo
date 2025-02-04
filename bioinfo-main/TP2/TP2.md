# TP2
## Lecture de fichiers
### 1.
On utilise la méthode read de SeqIO pour lire un fichier puis l'afficher avec seq
```
>>> from Bio import SeqIO
>>> record_genbank = SeqIO.read("./data/sequence.gb", "genbank")
>>> print(record_genbank.seq)
>>> record_fasta = SeqIO.read("./data/sequence.fasta", "fasta")
>>> print(record_fasta.seq)
```
### 2.
Il suffit d'utiliser la fonction str() pour faire passer en string la séquence. On peut ainsi vérifier que les 2 séquences sont les mêmes
```
>>> print(str(record_fasta.seq) == str(record_genbank.seq))
>>> True
```

### 3.
On utilise la méthode complement() de la classe Seq pour afficher la séquence complémentaire
```
print(str(record_fasta.seq.complement()))
```

### 4.
Il existe les attributs id, name, description ou encore seq. Les attributs sont communs aux entrées mais peuvent ne pas être remplis. Par exemple FASTA permet principalement de récupérer la séquence sans features alors que l'entrée GenBank permet de récupérer les features.

### 5.
La première feature est source. Pour récupérer sa position on peut utiliser l'attribut location de l'objet SeqFeature.
```
print(record_genbank.features[0].location)
```

