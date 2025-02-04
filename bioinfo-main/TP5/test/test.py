import os, shutil, sys, unittest
sys.path.append('..')
from src.utils import *
from src.pwn import *
from Bio import motifs

# Pour tester la classe :
# python3 -m unittest -v test.py


class TestScript(unittest.TestCase):

    maxDiff = None

    def getPssm(self):
        """
            Permet de récupérer une matrice pssm pour les tests
        """
        matrice_fpm = None
        message_erreur = "La matrice récupérée n'est pas une matrice PSSM"
        with open(self.chemin_matrices) as handle:
            for matrice in motifs.parse(handle, "jaspar"):
                if matrice.matrix_id[0:6] == self.nom_matrice or matrice.matrix_id == self.nom_matrice:
                    matrice_fpm = matrice
        self.assertIsNotNone(matrice_fpm)
        return pwn2pssm(matrice_fpm.counts, 0.1)

    def setUp(self):
        """
            Initialisation des attributs utilisés dans les tests
        """
        self.liste = ["NM_007389", "NM_001100", "NM_002469"]
        self.longueur_promoteur = 30
        self.repertoire = "./data_test"
        self.nom_matrice = "MA0114"
        self.longueur_matrice = 9
        self.chemin_matrices = "./pfm/pfm.txt"
        # Les séquences ne correspondent pas aux arn, c'est uniquement pour les tests
        self.sequence = "TTTGCCAAAGTCCACA"
        self.dict_seq = {
            self.liste[0] : "TCGGCTAATGTCCACACAAATTCCTGTTT",
            self.liste[1] : "TACGCCAAAGTCCACATAAATTCCAGTCA",
        }
    
    def test_mrna_to_gene(self):
        """
            Vérifie que l'id du gene récupéré à partir d'un
            numéro d'accession d'un ARN messager est bon
        """
        self.assertEqual(mrna_to_gene(self.liste[0]), "11435")

    def test_download_promotors(self):
        """
            Vérifie que les fichiers sont créés avec les séquences
        """
        if os.path.isdir(self.repertoire):
            shutil.rmtree(self.repertoire)
        os.mkdir(self.repertoire)
        download_promotors(self.liste, self.longueur_promoteur, self.repertoire)
        for arn in self.liste:
            nom_fichier = arn + "_" + str(self.longueur_promoteur) + ".fa"
            chemin = os.path.join(self.repertoire, nom_fichier)
            self.assertTrue(os.path.isfile(chemin))
    
    def test_pwn2pssm(self):
        """
            Vérifie qu'on récupère bien une matrice PSSM
            à partir de la matrice fpm
        """
        matrice_fpm = None
        message_erreur = "La matrice récupérée n'est pas une matrice PSSM"
        with open(self.chemin_matrices) as handle:
            for matrice in motifs.parse(handle, "jaspar"):
                if matrice.matrix_id[0:6] == self.nom_matrice or matrice.matrix_id == self.nom_matrice:
                    matrice_fpm = matrice
        self.assertIsNotNone(matrice_fpm)
        matrice_pssm = pwn2pssm(matrice_fpm.counts, 0.1)
        self.assertIsInstance(matrice_pssm, motifs.matrix.PositionSpecificScoringMatrix, message_erreur)
    
    def test_scan_sequence(self):
        """
            Vérifie que les résultats du scan de la séquence sont cohérents
        """
        matrice_pssm = self.getPssm()
        retour_scan = scan_sequence(matrice_pssm, self.sequence, 0)
        retour_scan[0] = (retour_scan[0][0], round(retour_scan[0][1], 2))
        self.assertEqual(str([(5, 16.42)]), str(retour_scan))
    
    def test_scan_all_sequences(self):
        """
            Vérifie que les résultats du scan de toutes
            les séquences sont cohérents
        """
        matrice_pssm = self.getPssm()
        retour_scan = scan_all_sequences(matrice_pssm, self.dict_seq, 0)
        for key, value in retour_scan.items():
            for i in range(len(retour_scan[key])):
                retour_scan[key][i] = (retour_scan[key][i][0], round(retour_scan[key][i][1], 2))
        self.assertEqual(str({
                self.liste[0]: [(5, 6.53), (16, 8.26)], 
                self.liste[1]: [(5, 16.42), (16, 7.27)]
            }), str(retour_scan))
    
    def test_score_window(self):
        """
            Vérifie que la fenêtre obtenue est cohérente
            avec les bons scores et positions
        """
        matrice_pssm = self.getPssm()
        dict_arnm_couple = scan_all_sequences(matrice_pssm, self.dict_seq, 0)
        window = score_window(dict_arnm_couple, 5, 20, self.longueur_matrice, self.longueur_promoteur)

        self.assertEqual(window, (11.47, -26, -11, 2, {
                self.liste[0] : [(-26, 6.53)],
                self.liste[1] : [(-26, 16.42)]
            }))
    
    def test_score_window_to_json(self):
        """
            Vérifie que les résultats du test précédent
            sont bien passés au format JSON
        """
        matrice_pssm = self.getPssm()
        dict_arnm_couple = scan_all_sequences(matrice_pssm, self.dict_seq, 0)
        window = score_window(dict_arnm_couple, 5, 20, self.longueur_matrice, self.longueur_promoteur)
        json = score_window_to_json(window)
        self.assertEqual(json, {
            "score": 11.47,
            "debut": -26,
            "fin": -11,
            "nb_occurence": 2,
            "occurences": [
                {
                    "arn": self.liste[0],
                    "position": -26,
                    "score": 6.53
                },
                {
                    "arn": self.liste[1],
                    "position": -26,
                    "score": 16.42
                }
            ]
        })
    
    def test_window_threshold(self):
        """
            Vérifie que les résultats obtenus en analysant un à un
            les fenêtres sont au bon format et ceux attendus
        """
        matrice_pssm = self.getPssm()
        dict_arnm_couple = scan_all_sequences(matrice_pssm, self.dict_seq, 0)
        window = window_threshold(dict_arnm_couple, self.longueur_matrice, self.longueur_promoteur, 21, 0)
        self.assertEqual(window, 
            {
                1: {
                    'score': 11.47,
                    'debut': -31,
                    'fin': -11,
                    'nb_occurence': 2,
                    'occurences': [{'arn': self.liste[0], 'position': -26, 'score': 6.53},
                                    {'arn': self.liste[1], 'position': -26, 'score': 16.42}]
                    },
                2: {
                    'score': 9.62,
                    'debut': -27,
                    'fin': -7,
                    'nb_occurence': 4,
                    'occurences': [{'arn': self.liste[0], 'position': -26, 'score': 6.53},
                                    {'arn': self.liste[0], 'position': -15, 'score': 8.26},
                                    {'arn': self.liste[1], 'position': -26, 'score': 16.42},
                                    {'arn': self.liste[1], 'position': -15, 'score': 7.27}]
                    },
                3: {
                    'score': 7.76,
                    'debut': -25,
                    'fin': -5,
                    'nb_occurence': 2,
                    'occurences': [{'arn': self.liste[0], 'position': -15, 'score': 8.26},
                                    {'arn': self.liste[1], 'position': -15, 'score': 7.27}]
                    }
            }
        )

    if __name__ == '__main__':
        unittest.main()
