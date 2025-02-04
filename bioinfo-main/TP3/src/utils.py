from Bio import SeqIO
from Bio import Entrez

def find_cds(record):
    ret = []
    for feature in record.features:
        if feature.type == "CDS":
            location = feature.location
            ret.append((int(location.start), int(location.end)))
    return ret

def mrna_to_gene(numero_accession):
    Entrez.email = "gaetan.oudshoorn.etu@univ-lille.fr"
    handle_elink_gb = Entrez.elink(dbfrom="nucleotide", db="gene", id=numero_accession, linkname="nucleotide_gene")
    record_elink_gb = Entrez.read(handle_elink_gb)
    handle_elink_gb.close()
    try:
        result = record_elink_gb[0]["LinkSetDb"][0]["Link"][0]["Id"]
        # result = record_elink_gb[0]["IdList"][0]
    except:
        raise ValueError
    return result