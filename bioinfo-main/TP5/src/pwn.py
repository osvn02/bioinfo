import numpy as np


def pwn2pssm(matrice_fpm, pseudo_poids):
    """
        Transforme une matrice fpm (FrequencyPositionMatrix) 
        en pwn (PositionWeightMatrix) puis en matrice pssm (PositionSpecificScoringMatrix)
    """
    matrice_pwm = matrice_fpm.normalize(pseudo_poids)
    pssm = matrice_pwm.log_odds()
    return pssm

def scan_sequence(matrice_pssm, sequence, seuil_score):
    """
        Scan toute une séquence et récupère tous les scores avec leur
        position dépassant le seuil requis à partir d'une matrice pssm
    """
    pos_scores = []
    longueur_pssm = len(matrice_pssm.consensus)
    longueur_sequence = len(sequence)
    for i in range(longueur_sequence - longueur_pssm + 1):
        seq = sequence[i : i+longueur_pssm]
        score = matrice_pssm.calculate(seq)
        if score >= seuil_score:
            pos_scores.append((i, score))
    return pos_scores

def scan_all_sequences(matrice_pssm, dict_arn_sequence, seuil_score):
    """
        Scan toutes les séquences dans un dictionnaire (nom_arnm, sequence) et récupère
        leur score dépassant le seuil requis à partir d'une matrice pssm.
        Le dictionnaire récupéré a pour clé l'ARN et valeur la liste position-score
    """
    dict_seq_couple = dict()
    for arn in dict_arn_sequence.keys():
        dict_seq_couple[arn] = scan_sequence(matrice_pssm, dict_arn_sequence[arn], seuil_score)
    return dict_seq_couple

def score_window(dict_arnm_couple, debut, fin, longueur_scan_sequence, taille_promoteur):
    """
        Calcul le score moyen d'une fenêtre donnée et renvoie
        le résultat de cette fenêtre avec:
        score = score moyen de la fenêtre calculé à partir de chaque séquence
        debut = position de début de la fenêtre (par rapport au début du gene)
        fin = position de fin de la fenêtre (par rapport au début du gene)
        compteur = nombre de séquence analysé dans la fenêtre
        dict_arn_score = dictionnaire associant chaque arnm à ses positions et scores trouvés

        Les paramètres utilisés sont :
        dict_arnm_couple = qui a en clé un numéro d'accession d'ARNm et en valeur une liste
        de position-score
        longueur_scan_sequence = longueur de la matrice et donc du scan de la matrice sur 
        un interval dans la fenêtre. 
        taille_promoteur = taille de la séquence promotrice
        debut et fin = position de debut et de fin de la fenêtre
    """
    dict_arn_score = {}
    moyenne_globale = 0
    compteur = 0
    # Parcours chaque ARNm
    for key in dict_arnm_couple.keys():
        # Parcours les positions-scores de chaque ARNm
        for pos_score in dict_arnm_couple[key]:
            # -2 car debut et fin comptent tous les 2 pour 1
            if pos_score[0] >= debut and pos_score[0] < fin - (longueur_scan_sequence - 2):
                # Initialise la liste position-score pour l'arn si pas encore fait
                if not key in dict_arn_score:
                    dict_arn_score[key] = []
                moyenne_globale += pos_score[1]
                compteur += 1

                # Ajoute l'occurence au dictionnaire
                new_position = pos_score[0] - taille_promoteur - 1
                score_arrondi = round(np.float64(pos_score[1]), 2) # Il y a un erreur si c'est un float32
                dict_arn_score[key].append((new_position, score_arrondi)) 
    res = None
    # Si au moins une occurence a été trouvée, on renvoie le score de la fenêtre
    if moyenne_globale != 0 and compteur > 0:
        moyenne_globale /= compteur
        debut = debut - taille_promoteur - 1
        fin = fin - taille_promoteur - 1
        moyenne_globale = round(moyenne_globale, 2)
        res = (moyenne_globale, debut, fin, compteur, dict_arn_score)        
    return res

def window_threshold(dict_arnm_couple, longueur_scan_sequence, taille_promoteur, taille_window, seuil_window, min_occ):
    """
        Renvoie toutes les fenêtres qui respectent le seuil de fenêtre demandé
        sous le format JSON (sans aucun doublon)
    """
    resultat = dict()
    compteur = 1
    # Parcours les promoteurs
    for i in range(taille_promoteur - longueur_scan_sequence):
        # Résultat de la fenêtre
        res = score_window(dict_arnm_couple, i, i + taille_window - 1, longueur_scan_sequence, taille_promoteur)
        # Vérifie si la fenêtre est valide (existe, seuil, minimum d'occurence)
        if res != None and res[0] >= seuil_window and res[3] >= min_occ:
            res = score_window_to_json(res)
            # Vérifie si les occurences sont les mêmes que la fenêtre précédente
            if compteur - 1 in resultat:
                if res["occurences"] != resultat[compteur - 1]["occurences"]:
                    resultat[compteur] = res
                    compteur += 1 
            else:
                resultat[compteur] = res
                compteur += 1 
                               
    return resultat

def score_window_to_json(score_window):
    """
        Utilise le résultat de la fonction score_window pour en créer un JSON 
    """
    occurences = []
    for arn, occurence in score_window[4].items():
        for pos, score in occurence:
            occurences.append({
                "arn": arn,
                "position": pos,
                "score": score
            })
    score_window = {
        "score": score_window[0],
        "debut": score_window[1],
        "fin": score_window[2],
        "nb_occurence": score_window[3],
        "occurences": occurences
    }
    return score_window

def print_windows(window_threshold_result):
    """
        Affichage de toutes les fenêtres 
    """
    nom_facteur = "Ahr::Arnt"
    num_window = 1
    while num_window in window_threshold_result:
        window = window_threshold_result[num_window]
        moyenne = str(window["score"])
        debut = str(window["debut"])
        fin = str(window["fin"])
        nb_occurence = str(window["nb_occurence"])
        string = str(num_window) + " " + nom_facteur + " [" + debut + ':' + fin + "] " + moyenne + " " + nb_occurence
        print(string)
        for occurence in window["occurences"]:
            print('\t', num_window, occurence["arn"], nom_facteur, occurence["position"], occurence["score"])
        num_window += 1
