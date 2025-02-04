Oudshoorn Gaëtan
Villegas Navarro Octavio Salvador
# Installation et lancement
Pour utiliser ce projet, il suffit de l'importer en local depuis ce gitlab et d'avoir python3 et biopython d'installer ```pip3 install```.  
Le script s'utilise avec la commande suivante qui peut être modifiée :
```
python3 putativeTFBS.py -m ./data/pfm.txt -a MA0114 -t 0 -l 1000 -w 40 -s 0 -p 0.1 -mo 5 -o ./export.json NM_001100 NM_002469 NM_002470 NM_003279 NM_003281
```

Cette commande contient 9 paramètres suivis des numéros d’accessions des ARN à analyser.

|  Paramètre         | Alias | Description | Obligatoire |
| :----------------- | :---: | :---------- | :---------: |
| --pfm              | -m    | Chemin relatif du fichier contenant les matrices | X |
| --matrice          | -a    | Nom de la matrice utilisée | X |
| --threshold        | -t    | Seuil de score | X |
| --promoter-length  | -l    | Longueur du promoteur | |
| --window-size      | -w    | Longueur de la fenêtre | |
| --window-threshold | -s    | Seuil de score de fenêtre | X |
| --pseudocount      | -p    | Pseudo-poids à ajouter à la matrice | |
| --min-occ          | -mo   | Minimum d'occurence par fenêtre | |
| --output           | -o    | Nom du fichier pour l'export en JSON | X |


Il est possible d'avoir les informations sur les paramètres avec la commande ```python3 putativeTFBS.py --help```

# Informations sur les dossiers
Dans ce projet il existe 3 dossiers :
- ./data/
- ./src/
- ./test/

## ./data/
Ce dossier contient toutes les données qui seront utilisées. C'est-à-dire le fichier pfm.txt, qui contient les matrices jaspar qui ont été récupérées sur https://jaspar.elixir.no/.
De plus ce dossier contiendra tous les promoters qui seront téléchargés par le script sous la forme [Numero_Accession_ARN]_[Longueur_Sequence].fa.

## ./src/
Ce dossier contient 2 fichiers **utils.py** et **pwn.py**.
#### utils.py
Ce fichier contient 2 fonction **mrna_to_gene** et **download_promotors**.  
Le premier permet de récupérer à partir d'un numéro d'accession d'un ARN messager l'id du gene correspondant. Cet id sera utilisé dans la 2ème fonction qui permet de télécharger des séquences promotrices d'une longueur et dans un répertoire qui sont passés en paramètre (par défaut le téléchargement se fait dans le repertoire courant du script). 

#### pwn.py
Ce fichier contient toutes les fonctions qui permettent de manipuler des matrices, de scanner des séquences et des fenêtres en passant par l'affichage des fenêtres pour le script. 

## ./test/
Ce dossier contient un fichier **test.py** qui contient une classe de test **TestScript** qui est la classe qui teste toutes les fonctions dans **pwn.py** et **utils.py**. De plus il y a 2 sous-dossiers, **./pfm/** qui contient uniquement le fichier **pfm.txt** (avec les matrices jaspar) et **./data_test/** qui est le dossier où le test de la fonction **download_promotors** télécharge les fichiers qui peuvent être vérifié.

# Présentation des résultats
## JSON
Le JSON récupéré avec toutes les fenêtres valides via script est de la forme suivante :
```json
{
    "1": {
        "score": 2.48,
        "debut": -240,
        "fin": -201,
        "nb_occurence": 5,
        "occurences": [
            {
                "arn": "NM_001100",
                "position": -231,
                "score": 0.91
            },
            {
                "arn": "NM_002469",
                "position": -209,
                "score": 4.42
            },
            {
                "arn": "NM_002470",
                "position": -225,
                "score": 1.78
            },
            {
                "arn": "NM_003279",
                "position": -218,
                "score": 2.56
            },
            {
                "arn": "NM_003281",
                "position": -223,
                "score": 2.73
            }
        ]
    }
}
```