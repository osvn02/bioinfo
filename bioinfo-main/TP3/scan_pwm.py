import argparse
from src.pwn import *
from src.utils import *
from Bio import SeqIO
from Bio import Entrez
from Bio import motifs

def main():
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fichier_jaspar', type=str, help='Chemin du fichier qui contient les matrices Jaspar')
    parser.add_argument('genbank_id', type=str, help='Identifiant Gankbank d\'un ARNm')
    parser.add_argument('taille_promoter', type=int, help='Taille de la séquence promotrice')
    parser.add_argument('seuil_score', type=float, help='Seuil de score')
    args = parser.parse_args()
    
    id_matrice = "MA0083"
    #id_matrice = "MA0114"
    #id_matrice = "MA0056"

    matrice_fpm = ''
    with open(args.fichier_jaspar) as handle:
        for matrice in motifs.parse(handle, "jaspar"):
            if matrice.matrix_id[0:6] == id_matrice or matrice.matrix_id == id_matrice:
                print(matrice.matrix_id, id_matrice)
                matrice_fpm = matrice
                print(matrice_fpm)
    
    matrice_pssm = pwn2pssm(matrice_fpm.counts, 0.1)
    print(matrice_pssm)
    # Récupère le gene
    gene_id = str(mrna_to_gene(args.genbank_id))
    print(gene_id)

    # Récupère le génome
    handle_gene_id = Entrez.esummary(db="gene", id=gene_id, retmode="text")
    record_gene_id = Entrez.read(handle_gene_id)
    handle_gene_id.close()
    genome_info = record_gene_id["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
    print("Info genome :")
    print(genome_info, '\n')
    chromosome_id = genome_info['ChrAccVer']
    chrstart_amont = int(genome_info["ChrStart"]) + 2
    # ChrStop : 73393624 (où 25 est le dernier)
    # ChrStart : 73410681 (où 82 est le premier)
    # Pour récupérer la séquence en amont il faut prendre
    # les chromosomes en commençant par le start + 2 (ici droite à gauche)
    # car on est sur le brin positif

    # Récupère la séquence en amont (ici on est sur le brin positif)
    handle_seq_amont = Entrez.efetch(db="nucleotide", id=chromosome_id, rettype="gb", retmode="text", seq_start=chrstart_amont, seq_stop=(chrstart_amont + args.taille_promoter - 1), strand=2)
    record_seq_amont = SeqIO.read(handle_seq_amont, 'gb')
    print("Séquence récupérée :\n")
    print(record_seq_amont.seq, '\n\n')
    scan = scan_sequence(matrice_pssm, record_seq_amont.seq, args.seuil_score)
    print(scan)
    return scan

if __name__ == "__main__":
    main()