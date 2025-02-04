import os
from Bio import SeqIO
from Bio import Entrez

def mrna_to_gene(numero_accession):
    """
        Récupère l'ID du gene à partir du numero d'accession d'un ARN.
    """
    Entrez.email = "gaetan.oudshoorn.etu@univ-lille.fr"
    handle_elink_gb = Entrez.elink(dbfrom="nucleotide", db="gene", id=numero_accession, linkname="nucleotide_gene")
    record_elink_gb = Entrez.read(handle_elink_gb)
    handle_elink_gb.close()
    try:
        result = record_elink_gb[0]["LinkSetDb"][0]["Link"][0]["Id"]
    except:
        raise ValueError
    return result

def download_promotors(list_arnm, len_seq_promo, repertoire="."):
    """
        Télécharge toutes les séquences promotrices depuis une liste d'ARN et 
        une longueur de séquence dans un repertoire spécifique.
        Les fichiers téléchargés sont sous la forme :
        [Numero_Accession_ARN]_[Longueur_Sequence].fa
    """
    for arn in list_arnm:
        # Récupère le gène
        gene_id = str(mrna_to_gene(arn))

        # Récupère le génome
        handle_gene_id = Entrez.esummary(db="gene", id=gene_id, retmode="text")
        record_gene_id = Entrez.read(handle_gene_id)
        handle_gene_id.close()
        genome_info = record_gene_id["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
        chromosome_id = genome_info['ChrAccVer']
        chrstart_amont = int(genome_info["ChrStart"]) + 2

        # Récupère la séquence en amont
        handle_seq_amont = Entrez.efetch(db="nucleotide", id=chromosome_id, rettype="fasta", retmode="text", seq_start=chrstart_amont, seq_stop=(chrstart_amont + len_seq_promo - 1), strand=2)
        record_seq_amont = SeqIO.read(handle_seq_amont, 'fasta')
        # Enregistrer les séquences
        nom_fichier = arn + "_" + str(len_seq_promo) + ".fa"
        chemin = os.path.join(repertoire, nom_fichier)
        with open(chemin, 'w') as fichier:
            fichier.write(str(record_seq_amont.seq))
    return
