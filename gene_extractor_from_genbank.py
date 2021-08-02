from Bio import SeqIO #This is useful for reading the input file
from Bio.Seq import Seq#This is useful for making the sequences in the end of the program 
from Bio.SeqRecord import SeqRecord#This is helpful for making the sequence records
###################################################################################
record=SeqIO.parse('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22.gbk','genbank')#input genebank file
for each in record:#each record is captured here
    feature=(each.features)#feature of each record is captured here
    gene_annotation=each.annotations['source']#we capture the name of the sequence here;if you want to capture anything alongwith the name remove the boxbrackets and then see your output
    nucleotides=each.seq#here we capture everything in the sequence(raw)
    for records in feature:#from this loop we capture the features to make the sequence
        gene_product=str(records.qualifiers.get('product')).strip("']['")#to know what the gene does we use this line
        gene_id=str(records.qualifiers.get('protein_id')).strip("']['")#to know the id of the protein it produces
        sequence=str(records.extract(nucleotides))#here we make the annotated sequences on the basis of features(refer line-7 &9)
        if (gene_product=='None'):#(gene_id=='None' and gene_product=='None'):#sometimes you get none in output,this line helps avoiding it
            continue
        else:
            sequence_record=SeqRecord(Seq(sequence),id=str(gene_id),description=(gene_product+' '+gene_annotation))#defining the whole identification of the gene
            with open('/home/skd/Desktop/MAN_LIPI/PS_22/PS_22_nuc.fasta','a') as file1:#saving it in a file
                file1.write(sequence_record.format('fasta'))#assigning the format of the file
##################################################################################