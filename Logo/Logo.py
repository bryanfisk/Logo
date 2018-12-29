import sys

threshold = None
input_file = None
#argv = ['test', '-t', '20', '-i', 'msa.txt']

def print_help():
    print('''
        -h, -help       print this help
        -i, -input      input file
        -t, -threshold  set minimum consensus threshold (percentage)
        ''')

#get input
exp_input = False
try:
    for i in range(1, len(sys.argv)):
        if exp_input == True:
            exp_input = False
            continue
        if sys.argv[i] in ('-h', '-help'):
            print_help()
            sys.exit(2)
        elif sys.argv[i] in ('-i', '-input'):
            exp_input = True
            input_file = sys.argv[i + 1]           
        elif sys.argv[i] in ('-t', '-threshold'):
            exp_input = True
            threshold = int(sys.argv[i + 1])
        else:
            print("Unhandled option: ", sys.argv[i])
            sys.exit()

    output_width = 80
    sequence_count = 0
    consensus_list = []
    position_dictionaries = {}

    #get input from file
    with open(input_file, 'r') as file:
        inputstring = file.read()
        inputstring = inputstring.strip()
        input = inputstring.split('\n')
        file.close()

    #counting how many sequences there are by the fasta headers
    for line in input:
        if line[0] == '>':
            sequence_count += 1
    alignment_list = [''] * sequence_count

    #checks first letter to see if it's a header, if not header then concatenate line in alignment file to alignment_list
    for line in input:
        if line[0] == '>':
            sequence_count -= 1
        else:
            alignment_list[sequence_count - 1] = alignment_list[sequence_count - 1] + line

    #making dictionaries for each position that holds how many occurances of each residue there are at that position
    for position in range(len(alignment_list[0])):
        temp_dict = {}
        for sequence in range(len(alignment_list)):
            if alignment_list[sequence][position] in temp_dict:
                temp_dict[alignment_list[sequence][position]] += 1
            else:
                temp_dict[alignment_list[sequence][position]] = 1
        position_dictionaries['position' + str(position + 1)] = temp_dict

    #converting occurances to percentages
    for position in position_dictionaries:
        for letter in position_dictionaries[position]:
            position_dictionaries[position][letter] = position_dictionaries[position][letter] / len(alignment_list) * 100

    #replacing below threshold positions and gaps with ellipses
    for position in position_dictionaries:
        maximum = max([list(position_dictionaries[position].items())[k][1] for k in range(len(list(position_dictionaries[position].items())))])
        for item in list(position_dictionaries[position].items()):        
            if item[1] >= threshold and item[0] != '-':
                consensus_list.append(item[0])
                break
            elif item[1] < threshold and item[0] != '-':
                consensus_list.append('…')
                break

    #making a list of the index of all multi-ellipsis stretches
    consolidate_elipsis = [item for item in range(0, len(consensus_list) - 1) if consensus_list[item] == consensus_list[item + 1] and consensus_list[item] == '…']

    #making the indexing number lists to number residue positions
    numbers_list = [[str(i).ljust(4)[0] for i in range(1, len(alignment_list[0]) + 1)], 
                    [str(i).ljust(4)[1] for i in range(1, len(alignment_list[0]) + 1)], 
                    [str(i).ljust(4)[2] for i in range(1, len(alignment_list[0]) + 1)],
                    [str(i).ljust(4)[3] for i in range(1, len(alignment_list[0]) + 1)]]

    #removing multi-ellipsis stretches
    for f, k in enumerate(consolidate_elipsis):
        consensus_list.pop(k - f)
        for i in numbers_list:
            i.pop(k - f)

    #removing numbers for the ellipses
    remove_numbers_under_elipsis = [item for item in range(len(consensus_list)) if consensus_list[item] == '…']
    for k in remove_numbers_under_elipsis:
        for i in numbers_list:
            i.pop(k)
            i.insert(k, ' ')

    #print out output
    for i in range((len(consensus_list[0]) // output_width) + 1):
        print(''.join(consensus_list[i * output_width : (i + 1) * output_width]))
        for numlist in numbers_list:
            print(''.join(numlist[i * output_width : (i + 1) * output_width]))

except SystemExit:
    pass
except TypeError:
    if input_file == None:
        print("Input file is required. (-i inputfile)")
    if threshold == None:
        print("Threshold is required. ( -t threshold)")

