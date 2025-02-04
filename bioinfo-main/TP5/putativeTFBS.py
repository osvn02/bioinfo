import json
import argparse
from src.utils import *
from src.pwn import *
from Bio import SeqIO
from Bio import Entrez
from Bio import motifs

# Commande utilisée :
# python3 putativeTFBS.py -m ./data/pfm.txt -a MA0114 -t 0 -l 1000 -w 40 -s 0 -p 0.1 -o ./export.json NM_001100 NM_002469 NM_002470 NM_003279 NM_003281

def main():
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--pfm', type=str, help='Chemin relatif du fichier contenant les matrices', required=True)
    parser.add_argument('-a', '--matrice', type=str, help='Nom de la matrice utilisée', required=True)
    parser.add_argument('-t', '--threshold', type=int, help='Seuil de score', required=True)
    parser.add_argument('-l', '--promoter-length', type=int, nargs='?', default=1000, help='Longueur du promoteur')
    parser.add_argument('-w', '--window-size', type=int, nargs='?', default=40, help='Longueur de la fenêtre')
    parser.add_argument('-s', '--window-threshold', type=int, help='Seuil de score de fenêtre', required=True)
    parser.add_argument('-p', '--pseudocount', type=float, nargs='?', default=0.1, help='Pseudo-poids')
    parser.add_argument('-mo', '--min-occ', type=int, nargs='?', default=0, help='Minimum d\'occurence par fenêtre')
    parser.add_argument('-o', '--output', type=str, help='Nom du fichier pour l\'export en JSON', required=True)
    parser.add_argument('ids_mrna', nargs='*', help='Les ARNm à télécharger et analyser')
    args = parser.parse_args()

    # Répertoire où seront stockés les séquences promotrices
    repertoire = "./data"

    # Télécharge tous les promoteurs des ANRm passé en argument
    download_promotors(args.ids_mrna, args.promoter_length, './data')

    # Crée un dictionnaire avec en clé le nom de l'ARNm et en valeur sa séquence
    dict_seq = dict()
    for arn in args.ids_mrna:
        nom_fichier = arn + "_" + str(args.promoter_length) + ".fa"
        chemin = os.path.join(repertoire, nom_fichier)
        with open(chemin, 'r') as fichier:
            dict_seq[arn] = fichier.read()

    # Récupère la matrice fpm à utiliser passé en argument
    matrice_fpm = ''
    with open(args.pfm) as handle:
        for matrice in motifs.parse(handle, "jaspar"):
            if matrice.matrix_id[0:6] == args.matrice or matrice.matrix_id == args.matrice:
                matrice_fpm = matrice

    # Créer la matrice PSSM à partir de la matrice fpm
    matrice_pssm = pwn2pssm(matrice_fpm.counts, args.pseudocount)
    longueur_scan_sequence = len(matrice_pssm[0])

    # Récupère un dictionnaire avec en clé le nom de l'ARNm
    # et en valeur une liste des couples position-score
    dict_arnm_couple = scan_all_sequences(matrice_pssm, dict_seq, args.threshold)

    # Récupération de toutes les fenêtres dont les scores
    # dépassent le seuil de score de fenêtre en récupérant
    # pour chaque fenêtre, les scores qui dépassent le seuil 
    window = window_threshold(dict_arnm_couple, longueur_scan_sequence, args.promoter_length, args.window_size, args.window_threshold, args.min_occ)
    
    # Affichage
    print_windows(window)

    # Créer l'objet JSON dont le fichier doit être créé
    json_object = json.dumps(window, indent=4)

    # Créer / écrit le JSON dans le fichier
    with open(args.output, "w") as fichier:
        fichier.write(json_object)

    return window

if __name__ == "__main__":
    main()