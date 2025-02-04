# TP3
## Calcul de score à partir de matrices de fréquences
### Q1
Pour trouver combien de matrices sont lus, on récupère avec la méthode parse de motifs toutes les matrices. On a juste à afficher la longueur de la liste retournée
```
with open("data/pfm.txt") as handle:
    retour = motifs.parse(handle, "jaspar")
print(len(retour))
```

### Q2
A partir d'une matrice FrequencyPositionMatrix, on peut utiliser la méthode normalize avec un pseudo-poids en paramètre pour créer une nouvelle matrice PositionWeightMatrix. 
```
matrice_corrigee = matrice.normalize(1)
```
où matrice est une matrice FrequencyPositionMatrix. 
On peut ajouter un pseudo_poids faible comme 0.1 pour ne pas avoir de valeurs à 0 tout en gardant quasiment la même proportion.

### Q3
Pour obtenir une matrice PSSM, on utilise la matrice avec le pseudo-poids créée et on utilise la méthode log_odds de la classe PositionWeightMatrix pour obtenir une matrice PositionSpecificScoringMatrix (PSSM).
```
matrice_pssm = matrice_corrigee.log_odds()
```

### Q4
Pour rechercher les occurences d'une PSSM dans une séquence on utilise la matrice PSSM et sa méthode calculate avec en paramètre l'objet Seq de la séquence
```
matrice_pssm.calculate(sequence_amont.seq)
```

## Mise en place de fonctions utiles pour manipuler les PWM
### Q1
Pour passer d'une matrice FrequencyPositionMatrix à une matrice PSSM, il suffit de passer d'abord par la matrice pseudo-poids pour ensuite créer la matrice PSSM comme précédemment
```
def pwn2pssm(matrice, pseudo_poids):
    matrice_pwm = matrice.normalize(pseudo_poids)
    pssm = matrice_pwm.log_odds()
    return pssm
```

### Q2
On parcourt la séquence de protéines et à chaque fois on récupère, un bout de la séquence. Le score PSSM de ce qui est récupéré est calculé en utilisant la méthode calculate de la matrice PSSM puis on vérifie que le score obtenu respecte le seuil pour l'ajouter à la liste des scores. Puis on renvoie la liste avec les positins et les motifs.

## Script pour le rendu intermédiaire
Pour ce script on entre en paramètre :  
- fichier_jaspar : le chemin vers le fichier jaspar qui contient les matrices
- genbank_id : l'identifiant genbank de l'ARNm
- taille_promoter : la taille de la séquence promotrice à récupérer
- seuil_score : le seuil de score

On récupère ensuite la matrice voulu dans le fichier jaspar avec le chemin en paramètre.
```
matrice_fpm = ''
with open(args.fichier_jaspar) as handle:
    for matrice in motifs.parse(handle, "jaspar"):
        if matrice.matrix_id[0:6] == id_matrice or matrice.matrix_id == id_matrice:
            print(matrice.matrix_id, id_matrice)
            matrice_fpm = matrice
            print(matrice_fpm)
```
Puis on utilise la fonction pwn2pssm pour faire passer la matrice en pssm avec le pseudo-poids.
```
matrice_pssm = pwn2pssm(matrice_fpm.counts, 0.1)
```
Ensuite on récupère le gene avec la fonction mrna_to_gene en utilisant en paramètre l'ID le l'ARNm. Ceci envoie l'ID link du gene puis on récupère les informations du gene avec esummary et read de Entrez.  
```
# Récupère le gene
gene_id = str(mrna_to_gene(args.genbank_id))
print(gene_id)
# Récupère le génome
handle_gene_id = Entrez.esummary(db="gene", id=gene_id, retmode="text")
record_gene_id = Entrez.read(handle_gene_id)
handle_gene_id.close()
```
Ainsi on peut ensuite récupérer les informations sur le génome dont le début, la fin et l'id du chromosome pour faire une requête efetch sur la base nucleotide. Dans cette requête on y met l'id du chromosome les indices de début et de fin de la séquence voulue et le brin que l'on veut récupérer (positif ou négatif).  
```
genome_info = record_gene_id["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
chromosome_id = genome_info['ChrAccVer']
chrstart_amont = int(genome_info["ChrStart"]) + 2
handle_seq_amont = Entrez.efetch(db="nucleotide", id=chromosome_id, rettype="gb", retmode="text", seq_start=chrstart_amont, seq_stop=(chrstart_amont + args.taille_promoter - 1), strand=2)
record_seq_amont = SeqIO.read(handle_seq_amont, 'gb')
```
Une fois qu'on a récupéré la séquence il suffit de la lire avec read et d'utiliser la fonction scan_sequence en passant en paramètre la matrice PSSM, la séquence en amont et le seuil de score. 
```
print(record_seq_amont.seq, '\n\n')
scan = scan_sequence(matrice_pssm, record_seq_amont.seq, args.seuil_score)
print(scan)
return scan
```
