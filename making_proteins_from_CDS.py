#Making the protein from CDS
import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import Bio.Data.CodonTable
from Bio.SeqRecord import SeqRecord
###############################################################################
def purifying_sequences():
    cds=os.listdir('/home/skd/Desktop/MAN_LIPI/PS_22')#input directory of the nucleotides
    os.mkdir('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds')#enter the the directory for all the output cds files
    for i in range(len(cds)):
        sequence=[]
        for each in SeqIO.parse('/home/skd/Desktop/MAN_LIPI/PS_22/'+cds[i],'fasta'):#input directory of the nucleotides
            if(len(each.seq)%3==0):
                sequence.append(SeqRecord(Seq(str(each.seq)),id=each.id,description=each.description))
        SeqIO.write(sequence,'/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds/'+cds[i],'fasta')#enter the the directory for the cds files generated in the multiples of 3
###############################################################################
def translating_sequences():
    cds=os.listdir('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds')#enter the the directory for the cds files generated in the multiples of 3
    os.mkdir('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_protein')#directory where you want your protein files
    for i in range(len(cds)):
        sequence=[]
        for each in SeqIO.parse('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds/'+cds[i],'fasta'):#enter the the directory for the cds files generated in the multiples of 3
            try:
            sequence.append(SeqRecord(Seq(str(((each.seq)).translate(table=11,cds=True,to_stop=True))),id=each.id,description=each.description))#from gbk file you can check for the codon table recommended and change the table=11 to other numbers
            except Bio.Data.CodonTable.TranslationError:
                continue
        SeqIO.write(sequence,'/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds/'+cds[i],'fasta')#enter the the directory for the cds files generated in the multiples of 3
###############################################################################
if __name__=='__main__':
    try:
        shutil.rmtree('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_all_%3==0_cds')#input directory of the nucleotides
        shutil.rmtree('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_protein')
    except OSError as e:
        print("Maybe I did it,but do check.Yo! bitcracker")
    purifying_sequences()
    translating_sequences()
################################################################################
