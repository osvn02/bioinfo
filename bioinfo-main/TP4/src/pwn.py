def pwn2pssm(matrice, pseudo_poids):
    matrice_pwm = matrice.normalize(pseudo_poids)
    pssm = matrice_pwm.log_odds()
    return pssm

def scan_sequence(matrice_pssm, sequence, seuil_score):
    # Garder brin positif
    scores = []
    longueur_pssm = len(matrice_pssm.consensus)
    longueur_sequence = len(sequence)
    
    for i in range(longueur_sequence - longueur_pssm + 1):
        seq = sequence[i : i+longueur_pssm]
        score = matrice_pssm.calculate(seq)
        if score >= seuil_score:
            scores.append((i, score))
    return scores

def scan_all_sequences(matrice_pssm, dict_sequence, seuil_score):
    dict_seq_couple = dict()
    for arn in dict_sequence.keys():
        dict_seq_couple[arn] = scan_sequence(matrice_pssm, dict_sequence[arn], seuil_score)
    return dict_seq_couple

def score_window(dict_arnm_couple, debut, fin, longueur_scan_sequence):
    dict_moyenne_par_arn = {}
    for key in dict_arnm_couple.keys():
        moyenne = 0
        compteur = 0
        for pos_score in dict_arnm_couple[key]:
            # -2 car debut et fin comptent tous les 2 pour 1
            if pos_score[0] >= debut and pos_score[0] < fin - (longueur_scan_sequence - 2):
                moyenne += pos_score[1]
                compteur += 1
                # print(key, pos_score)
        if moyenne != 0:
            moyenne /= compteur
            dict_moyenne_par_arn[key] = moyenne
    # print("\nMoyenne pour chaque ARNm :\n", dict_moyenne_par_arn)
    compteur = 0
    total = 0
    for key in dict_moyenne_par_arn.keys():
        total += dict_moyenne_par_arn[key]
        compteur += 1
    res = None
    if(compteur > 0):
        # print("\nMoyenne des ARNm :", total / compteur)
        res = (total / compteur, debut, fin)
    return res

def best_window(dict_arnm_couple, longueur_scan_sequence, taille_promoteur, taille_window):
    meilleur_resultat = None
    for i in range(taille_promoteur - longueur_scan_sequence):
        if meilleur_resultat != None:
            resultat = score_window(dict_arnm_couple, i, i + taille_window, longueur_scan_sequence)
            if resultat != None and resultat[0] > meilleur_resultat[0]:
                meilleur_resultat = resultat
        else:
            meilleur_resultat = score_window(dict_arnm_couple, i, i + taille_window, longueur_scan_sequence)
    return meilleur_resultat