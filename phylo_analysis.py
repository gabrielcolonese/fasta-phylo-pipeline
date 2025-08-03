from Bio.Seq import Seq, UndefinedSequenceError
from Bio import Blast
from Bio import AlignIO, SeqIO
from Bio import Phylo
from matplotlib import pyplot as plt
import sys
import re
import xml.etree.ElementTree as ET
import pandas as pd
import numpy as np

num_species = int(sys.argv[2])

def get_ORFs(seq):
    orf_list = []
    reverse = seq.reverse_complement()
    for i in range(0, 3):
        orf_list.append(seq[i:].translate())
        orf_list.append(reverse[i:].translate())
    return orf_list

    
def screen1(fname):
    print("Type 0 if the file is a protein sequence, 1 if it is a DNA sequence:")
    opt = int(input())

    for record in SeqIO.parse(fname, "fasta"):
        my_seq = record.seq

    write_to_report("Input sequence", text=f"The input sequence was {my_seq}")

    #if DNA, get ORFs and ask the user to select the one they want:
    if opt:
        aa_seq = screen_orf(my_seq)
    else:
        aa_seq = my_seq
    
    return blast_search(aa_seq)


def screen_orf(seq):
    aa_seq_list = get_ORFs(seq)
    print("Type the number that corresponds to the intended AA sequence:")
    j = 1
    for i in range(0, len(aa_seq_list), 2):
        print(f'{i+1} (Forward Strand {j}): {aa_seq_list[i]}\n')
        print(f'{i+2} (Reverse Strand {j}): {aa_seq_list[i]}\n')
        j += 1
    opt = int(input())-1

    #make error exception (try-catch) if value is not in the 1:6 range
    
    # return blast_search(aa_seq_list[opt])
    output = f"The selected ORF was {aa_seq_list[opt]}"
    write_to_report("Selected ORF", output)
    return aa_seq_list[opt]


def blast_search(my_seq):
    print("Commencing BLAST...\n")
    #result = Blast.qblast("blastp", database="nr", sequence=my_seq, format_type="Text")
    result = Blast.qblast("blastp", database="nr", sequence=my_seq, format_type="XML")
    with open("my_blast.xml", "wb") as out_stream:
        out_stream.write(result.read())
    result.close()
    print("BLAST done!!")

    blast_record = read_blast("my_blast.xml")
    return write_top_k(blast_record, num_species)


def read_blast(fname):
    # result_stream = open("my_blast.xml", "rb")
    result_stream = open(fname, "rb")
    #blast_records = Blast.parse(result_stream)
    blast_records = Blast.read(result_stream)
    #result_stream.close()
    return blast_records


#Writes top k files (sorted by e-value) in a fasta file
def write_top_k(blast_record, k):
    names = []
    print("Preparing to write the file for the MSA\n")
    arq = open("my_MSA_in.fasta", "w")

    #sorts the blast output based on score
    sort_key = lambda hit: hit[0].score
    blast_record.sort(key=sort_key, reverse=True)

    id_list = [] #list to be called when converting the xml files to pandas df
    list = []
    #k = len(blast_record)
    for i in range(0, len(blast_record)):
        if len(list) > k:
            break
        # if k > len(blast_record):
        #     print("Error: number of species inserted is larger than BLAST's output")
        #     sys.exit(0)
        #     break

        hit = blast_record[i]
        alignment = hit[0]

        regex = re.findall(r"(\[.+?\])+", alignment.target.description)
        for item in regex:
            splitted = item.split("[")[1]
            splitted = splitted.split("]")[0]
            splitted = splitted.split()
            item = splitted[0] + '_' + splitted[1]
            
            if item not in list and item != "synthetic_construct" and item != "Human_ORFeome":
                is_undefined = False
                try:
                    bytes(alignment.target.seq)
                except UndefinedSequenceError:
                    is_undefined = True

                if is_undefined:
                    continue

                # print(f"appending {item}")
                list.append(item)

                # print(f"id: {alignment.target.id}")
                # print(f"description: {alignment.target.description}")
                #fasta_header = ">" + alignment.target.id + alignment.target.description
                id_list.append(alignment.target.id + alignment.target.description)

                re_split = re.search(r"\|", alignment.target.id)
                id = alignment.target.id[re_split.start()+1 :]
                fasta_header = ">" + id + item
                # fasta_header = ">" + alignment.target.id + item
                fasta_header = fasta_header.replace("|", "_")
                fasta_seq = str(alignment.target.seq)
                arq.write(fasta_header + "\n")
                arq.write(fasta_seq + "\n")
                arq.write("\n")
                names.append(item)

                break



        #if i == 10: #for some reason i=10 throws a an UndefinedSequenceError
        ###if alignment.target.id == "gb|PNI99818.1|": #for some reason this one (i=10) throws a an UndefinedSequenceError
        ###    continue
        #print("\n")
        #print(alignment.target)
        #print(alignment.target.seq)
        # arq.write("\n")


    if len(names) == 0:
        sys.exit("Error: No relevant proteins found in the BLAST file")
    print("Done!")
    arq.close()

    write_to_report("Proteins selected out of the BLAST search", list = id_list)

    # print(names)
    return screen_MSA()


def screen_MSA():
    import subprocess
    print("Commencing MUSCLE")
    cmd = r'muscle-win64.v5.3 -align my_MSA_in.fasta -output my_MSA_out.txt'
    results = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
    #align = AlignIO.read("my_MSA_out.txt", "fasta")
    #print(align)
    return make_tree("my_MSA_out.txt")


def matrix_to_dataframe(mat):
    size = len(mat.names)
    # Create an empty matrix filled with NaNs or zeros
    full_matrix = np.full((size, size), np.nan)  # or use 0.0 instead of np.nan

    # Fill in the lower triangular part
    for i in range(size):
        full_matrix[i, :i+1] = mat.matrix[i]

    # Optionally mirror to the upper triangle if symmetric
    for i in range(size):
        for j in range(i+1, size):
            full_matrix[i][j] = full_matrix[j][i]

    # Create DataFrame
    df = pd.DataFrame(full_matrix, index=mat.names, columns=mat.names)
    return df


def make_tree(filename):
    msa = AlignIO.read(filename, "fasta")
    #names = [s.id for s in msa]
    #print(names)

    from Bio.Phylo import TreeConstruction
    calculator = TreeConstruction.DistanceCalculator(model='blosum62')
    distance_matrix = calculator.get_distance(msa) #get distance matrix
    df = matrix_to_dataframe(distance_matrix)
    write_to_report("Distance matrix", table=df)
    #write_to_report("Distance matrix", str(distance_matrix))
    print(type(distance_matrix))

    constructor = TreeConstruction.DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)
    return draw_tree(tree)


def draw_tree(tree):
    # print(str(tree))
    # print(type(tree))
    fig  = plt.figure(figsize=(30, 15))
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig("fig.jpg")
    write_to_report("Phylogenetic tree", img="fig.jpg")
    return 1


def write_to_report(header, text="", img="", table = None, list = None):
    with open("report.md", "a") as f:
        f.write("\n")
        f.write(f"# {header}\n\n")

        if text != "":
            f.write(f"{text}\n\n")
        if img != "":
            f.write(f'![Tree]({img})\n\n')
        if table is not None:
            f.write(table.to_markdown(tablefmt="github"))
        if list is not None:
            i = 1
            for item in list:
                f.write(f'{i}. {item}\n\n')
                i += 1
            
        

#User says if its DNA or protein seq
#if DNA show corresponding protein seqs for reading frames, if protein proceed to BLAST search
#database options available at https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome
# nr, swissprot, pdb, etc etc 


if __name__ == "__main__":
    filename = sys.argv[1]
    with open("report.md", "w") as f:
        f.write("")
    screen1(filename)

    # blast_record = read_blast("my_blast.xml")
    # write_top_k(blast_record, 10)