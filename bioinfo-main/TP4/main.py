import argparse
from src.utils import *
from src.pwn import *
from Bio import SeqIO
from Bio import Entrez
from Bio import motifs

# Télécharge tous les promoteurs des ANRm
def main():
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--threshold', type=int, help='Seuil de score')
    parser.add_argument('--promoter-length', type=int, nargs='?', default=1000, help='Longueur du promoteur')
    parser.add_argument('--window-size', type=int, nargs='?', default=40, help='Longueur de la fenêtre')
    parser.add_argument('ids_mrna', nargs='*', help='Les ARNm')
    args = parser.parse_args()
    print("Seuil de score :", args.threshold)
    print("Longueur du promoteur :", args.promoter_length)
    print("Longueur de la fenêtre :", args.window_size)
    print("Les ARNm :", args.ids_mrna)
    repertoire = "./data"

    download_promotors(args.ids_mrna, args.promoter_length, './data')

    # Crée un dictionnaire avec en clé le nom de l'ARNm
    # et en clé la séquence
    dict_seq = dict()
    for arn in args.ids_mrna:
        nom_fichier = arn + "_" + str(args.promoter_length) + ".fa"
        chemin = os.path.join(repertoire, nom_fichier)
        with open(chemin, 'r') as fichier:
            dict_seq[arn] = fichier.read()
    # Récupère la matrice
    id_matrice = "MA0114"
    matrice_fpm = ''
    with open("./data/pfm.txt") as handle:
        for matrice in motifs.parse(handle, "jaspar"):
            if matrice.matrix_id[0:6] == id_matrice or matrice.matrix_id == id_matrice:
                matrice_fpm = matrice
    matrice_pssm = pwn2pssm(matrice_fpm.counts, 0.1)
    longueur_scan_sequence = len(matrice_pssm[0])
    # Récupère un dictionnaire avec en clé le nom de l'ARNm
    # et en clé une liste des couples position-score
    dict_arnm_couple = scan_all_sequences(matrice_pssm, dict_seq, args.threshold)

    # Test de récupération de la meilleure fenêtre
    return best_window(dict_arnm_couple, longueur_scan_sequence, args.promoter_length, args.window_size)

if __name__ == "__main__":
    main()