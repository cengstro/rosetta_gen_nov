# sorts all fasta entries in `data/seq/genbank` and 'data/seq/single_cell` 
# which are organized according to Genewiz sequencing batch (along with their respective tracefiles)
# and sorts them by genetic marker: one of ITS2, rbcL, or 18S

from datetime import date
import glob
from Bio import SeqIO
import re

today = date.today().strftime("%Y-%m-%d")

# list of filepaths
my_fastas = glob.glob("/home/cengstro/ownCloud/proj/single_cell/data/seq/single_cell/[CE|BR|JA]*/*.fasta")
gb_fastas = glob.glob("/home/cengstro/ownCloud/proj/single_cell/data/seq/genbank/all_markers_for_tree.fasta")
all_fastas = my_fastas + gb_fastas
# print(all_fastas)
# print(len(all_fastas), len(my_fastas))

# name the three output files to be created (e.g. rbcL_Sep20.fasta)
its2_out_path = "/home/cengstro/ownCloud/proj/single_cell/data/seq/its2_" + today +".fasta"
x18s_out_path = "/home/cengstro/ownCloud/proj/single_cell/data/seq/x18s_" + today +".fasta"
rbcl_out_path = "/home/cengstro/ownCloud/proj/single_cell/data/seq/rbcL_" + today +".fasta"

its2_out_handle=open(its2_out_path,'w')
x18s_out_handle=open(x18s_out_path,'w')
rbcl_out_handle=open(rbcl_out_path,'w')

# loop through the fastas and sort the entries into the desired output files
for fasta in all_fastas:
    parsed_fasta = SeqIO.parse(open(fasta),'fasta')
    print("opening: " + fasta)
    for record in parsed_fasta:
        my_seq = str(record.seq)
        my_header = str(record.description) # id just gives the first space delim in the header
        # print(my_header)
        if re.search('ITS2|internal', my_header, re.IGNORECASE): 
            print("writing ITS2: " + my_header)
            its2_out_handle.write(">" + my_header + "\n") # write header
            its2_out_handle.write(my_seq + "\n") # write seq
        if re.search('18S', my_header, re.IGNORECASE): # don't use elif, since a single header might contain both its2 and 18s
            print("writing 18S: " + my_header, re.IGNORECASE)
            x18s_out_handle.write(">" + my_header + "\n") # write header
            x18s_out_handle.write(my_seq + "\n") # write seq
        if re.search('rbcL|ribulose', my_header, re.IGNORECASE):
            print("writing rbcL: " + my_header)
            rbcl_out_handle.write(">" + my_header + "\n") # write header
            rbcl_out_handle.write(my_seq + "\n") # write seq
