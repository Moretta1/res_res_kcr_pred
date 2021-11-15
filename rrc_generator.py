from subscript import *
from math import isnan

'''
input: residue-residue contact prediction result
each column meaning:
1st number = residue 1
2nd number = residue 2
3rd and 4th number = d1 and d2 from CASP format (legacy format requirement. not relevant to predictions here)
5th number = your prediction! A probability (between 0 and 1) that the residues are in contact. With 0 = no contact and 1 = contact

output: .fasta file of residue-residue contact peptide and residue-residue composition(RRC)
'''




# read the original sequence for contact extraction
originals_sequence = readFasta('residues.fasta')
originals_sequence = np.array(originals_sequence)

# mapping ids:
mapped = pd.read_csv('id_mapping.csv', header=0)
mapped = np.asarray(mapped)
# find the corresponding sequence of content file and original file:
win_size = 10
file = "MP007031_results"

residue, marked_line = findMarkedLine(file, win_size)
# here the marked_line might be a table, as it may contains several sites

# read the predicted results:
path = file + "/cont.txt"

content = pd.read_table(path, header=0)
content = np.array(content)[1:, :]
segments = list()
empty_list = list()
for i in range(0, len(marked_line)):
    lower = marked_line[i, 3]
    upper = marked_line[i, 4]
    seg = marked_line[i, 0]
    if type(seg) is not str and isnan(seg):
        print('nan')
        continue
    position = lower + 1 + win_size
    ty = marked_line[i, 1]
    added = list()
    print(i)

    for j in range(0, len(content)):
        arr = str.split(str(content[j, 0]), sep=' ')
        # find residue 1 and residue 2
        index1 = int(arr[0])-1
        index2 = int(arr[1])-1

        # find predicted possibility
        prob = float(str.split(arr[4], sep=' ')[0])

        # threshold here: 0.8
        if prob >= 0.8:
            if index1 <= lower <= index2 <= upper:
                added.append(index1)

            elif lower < index1 <= upper <= index2:
                added.append(index2)

    # remove redundant contact residues
    added_unique = list(set(added))

    if len(added_unique) != 0:
        empty_list.append(ty)
        for k in added_unique:
            seg = seg + str(residue[int(k)])
        empty_list.append(seg)
        print(len(empty_list))

result = np.asarray(empty_list)

# save the result as fasta file
output = "output.fasta"
with open(output, 'a') as f:
    for each in result:
        if '-' in each:
            continue

        label = '>' + str(each) if len(each) == 1 else each
        f.write(str(label))
        f.write('\n')

# read the output.fasta then encoding it as RRC feature
data = readFasta("output.fasta")

aac = AAC(data)

savetsv(aac, "/Users/apple/Documents/residue-residue-contact/aac.tsv")
