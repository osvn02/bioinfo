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