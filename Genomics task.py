from collections import defaultdict
import numpy
def read_single_txt(file_name):
    file = open(file_name, "r")
    sequence_length = int(file.readline())
    sequences = [sequence.rstrip() for sequence in file.readlines()]
    print("The sequence length is :")
    print(sequence_length)
    print("The sequences are :")
    print(sequences)
    return sequence_length, sequences
def read_single_fasta(file_name):
    file = open(file_name, "r")
    sequence_length = int(file.readline())
    sequences = []
    while True:
        line = file.readline()
        if line == "":
            break
        sequences.append(file.readline().rstrip())
        file.readline()
        file.readline()
    print("The sequence length is :")
    print(sequence_length)
    print("The sequences are :")
    print(sequences)
    return sequence_length, sequences
def create_de_bruijn_single(sequence_length, sequences):
    k_minus_one = sequence_length - 1
    de_bruijn = defaultdict(list)
    for sequence in sequences:
        de_bruijn[sequence[0:k_minus_one]].append(sequence[1:sequence_length])
    print(de_bruijn)
    return de_bruijn

def find_path_single(de_bruijn):
    de_bruijn_keys = set(de_bruijn.keys())
    de_bruijn_values = set(numpy.concatenate(list(de_bruijn.values())).flat)
    start = de_bruijn_keys.difference(de_bruijn_values).pop()
    end = de_bruijn_values.difference(de_bruijn_keys).pop()
    path = []
    change_key = start
    path.append(change_key)
    while change_key != end:
        if len(de_bruijn[change_key]) != 0:
            path.append(de_bruijn[change_key][0])
            change_key = de_bruijn[change_key].pop(0)
    print(path)
    return path
def assembly_single(path):
    assembly = ""
    assembly += path[0]
    del path[0]
    last_bases = [sequence[-1] for sequence in path]
    assembly += "".join(last_bases)
    print(assembly)
    return assembly

def read_paired_txt(file_name):
    file = open(file_name, "r")
    first_line_split = file.readline().split(' ')
    sequence_length = int(first_line_split[0])
    gap_length = int(first_line_split[1])
    paired_sequences = [sequence for sequence in file.readlines()]
    initial_sequences = [pair.split('|')[0] for pair in paired_sequences]
    terminal_sequences = [pair.split('|')[1].rstrip() for pair in paired_sequences]
    print("The sequence length is :")
    print(sequence_length)
    print("The gap length is :")
    print(gap_length)
    print("The initial sequences are :")
    print(initial_sequences)
    print("The terminal sequences are :")
    print(terminal_sequences)
    return sequence_length, gap_length, initial_sequences, terminal_sequences
def create_de_bruijn_paired(sequence_length,initial_sequences, terminal_sequences):
    k_minus_one = sequence_length - 1
    de_bruijn = defaultdict(list)
    for i in range(len(initial_sequences)):
        initial_k_mer = (initial_sequences[i][0:k_minus_one], terminal_sequences[i][0:k_minus_one])
        terminal_k_mer = (initial_sequences[i][1:sequence_length], terminal_sequences[i][1:sequence_length])
        de_bruijn[initial_k_mer].append(terminal_k_mer)
    print(de_bruijn)
    return de_bruijn
def find_path_paired(de_bruijn):
    de_bruijn_keys = set(de_bruijn.keys())
    de_bruijn_values = set([tuple(numpy.concatenate(item)) for item in de_bruijn.values()])
    start = de_bruijn_keys.difference(de_bruijn_values).pop()
    end = de_bruijn_values.difference(de_bruijn_keys).pop()
    path = []
    change_key = start
    path.append(change_key)
    while (change_key[0] != end[0]) & (change_key[1] != end[1]):
        if len(de_bruijn[change_key]) != 0:
            path.append(de_bruijn[change_key][0])
            change_key = de_bruijn[change_key].pop(0)
    print(path)
    return path

def assembly_paired(sequence_length, gap_length, path):
    first_bases_initial = [sequence[0][0] for sequence in path]
    first_bases_terminal = [sequence[1][0] for sequence in path]
    del first_bases_initial[-1]
    del first_bases_terminal[-1]
    prefix_string = "".join(first_bases_initial) + path[-1][0]
    suffix_string = "".join(first_bases_terminal) + path[-1][1]
    assembly = prefix_string + suffix_string[len(suffix_string) - (sequence_length + gap_length):]
    print(assembly)
    return assembly
def check_testcase_single(assembly):
    file = open('SingleReadOutput.txt','r')
    answer = file.readline().rstrip()
    if answer == assembly:
        return True
    else:
        return False
def check_testcase_paired(assembly):
    file = open('ReadPairsOutput.txt','r')
    answer = file.readline().rstrip()
    if answer == assembly:
        return True
    else:
        return False

type_of_paires=(input ("What type of pairing you want to insert\n{single or paired read?} \n"))
if(type_of_paires == "single"):
 #first creat graph for single_paires (((-------------------------*****----------------------------------------)))
 type_of_files = (input("What type of file you want to read\n{txt or fastq?} \n"))
 if (type_of_files == "txt"):
    sequence_length, sequences = read_single_txt('SingleReadsInput.txt')
 elif(type_of_files == 'fastq'):
     sequence_length, sequences = read_single_fasta('SingleReads.fastq')
 print("The created graph for single_reads sequence you entered is : ")
 de_bruijn = create_de_bruijn_single(sequence_length, sequences)
 print("The chosen eulerian walk path through the graph is : ")
 path = find_path_single(de_bruijn)
 print("The assembled sequence is : ")
 assembly = assembly_single(path)
 print(len(assembly))

                             #========================================================#

elif(type_of_paires == "paired"):
 # second creat graph for paired_paires (((-------------------------*****----------------------------------------)))
 sequence_length, gap_length, initial_sequences, terminal_sequences = read_paired_txt("ReadPairsInput.txt")
 print("The created graph for paired_reads sequence you entered is : ")
 de_bruijn = create_de_bruijn_paired(sequence_length,initial_sequences, terminal_sequences)
 print("The chosen eulerian walk path through the graph is :")
 path = find_path_paired(de_bruijn)
 print("The assembled sequence is : ")
 assembly = assembly_paired(sequence_length, gap_length, path)

else:

        print("written mistake!")
