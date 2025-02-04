from Bio import SeqIO
from Bio import Entrez
from Bio import motifs
from src.utils import find_cds, mrna_to_gene

# record_genbank = SeqIO.read("./data/sequence.gb", "genbank")
# print(record_genbank.seq)

# record_fasta = SeqIO.read("./data/sequence.fasta", "fasta")
# print(record_fasta.seq)
# print(str(record_fasta.seq) == str(record_genbank.seq))

# print(str(record_fasta.seq.complement()))
# print(record_fasta.id, record_genbank.id)
# print(record_fasta.name, record_genbank.name)
# print(record_fasta.description, record_genbank.description)

# print(record_genbank.features[3])
#for feature in record_genbank.features:
#    print("\n\n")
#    print(feature)

# print(find_cds(record_genbank))

Entrez.email = "gaetan.oudshoorn.etu@univ-lille.fr"
# handle_gb = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="gb", retmode="text")
# read_gb = SeqIO.read(handle_gb, "genbank")
# print(read_gb, '\n')

# handle_fasta = Entrez.efetch(db="nucleotide", id="NM_007389", rettype="fasta", retmode="text")
# read_fasta = SeqIO.read(handle_fasta, "fasta")
# record_genbank = SeqIO.read("./data/sequence.gb", "genbank")
# print(read_fasta)
# print(str(read_gb.seq) == str(read_fasta.seq))
# print(str(find_cds(read_gb)) == str(find_cds(record_genbank)))

# handle_elink_gb = Entrez.elink(dbfrom="nucleotide", db="gene", id="NM_007389", linkname="nucleotide_gene")
# record_elink_gb = Entrez.read(handle_elink_gb)
# handle_elink_gb.close()
# print(record_elink_gb[0]["LinkSetDb"][0]["Link"][0]["Id"])
# print(record_elink_gb[0]["IdList"][0])
# print(record_elink_gb)


gene_id = str(mrna_to_gene("NM_007389"))
print("Gene id :",gene_id)
handle_gene_id = Entrez.esummary(db="gene", id=gene_id, retmode="text")
record_gene_id = Entrez.read(handle_gene_id)
handle_gene_id.close()
genome_info = record_gene_id["DocumentSummarySet"]["DocumentSummary"][0]["GenomicInfo"][0]
print("Info genome :", genome_info, '\n')
chromosome_id = genome_info['ChrAccVer']
chrstart_amont = int(genome_info["ChrStart"]) + 2
# ChrStop : 73393624 (où 25 est le dernier)
# ChrStart : 73410681 (où 82 est le premier)
# Pour récupérer la séquence en amont il faut prendre
# les chromosomes en commençant par le start + 2 (ici droite à gauche)
handle_seq_amont = Entrez.efetch(db="nucleotide", id=chromosome_id, rettype="fasta", retmode="text", seq_start=chrstart_amont, seq_stop=(chrstart_amont + 999), strand=2)
record_seq_amont = SeqIO.read(handle_seq_amont, 'fasta')
record_fasta_amont = SeqIO.read("./data/seq_amont.fa", "fasta")
print(record_seq_amont.seq, '\n')
print(record_fasta_amont.seq)
print(str(record_fasta_amont.seq) == str(record_seq_amont.seq))