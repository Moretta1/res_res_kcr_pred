import sys
import pandas as pd
import re
import numpy as np
import math
from collections import Counter
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, auc
import os

def readFasta(file):
    if not os.path.exists(file):
        print('Error: "' + file + '" does not exist.')
        sys.exit(1)

    with open(file) as f:
        records = f.read()

    if re.search('>', records) is None:
        print('The input file seems not in fasta format.')
        sys.exit(1)

    records = records.split('>')[1:]
    myFasta = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ARNDCQEGHILKMFPSTWYV-]', '-', ''.join(array[1:]).upper())
        myFasta.append([name, sequence])
    return myFasta

def read_svm(file):
    encodings = []
    labels = []
    with open(file, encoding='utf-16') as f:
        # with open(file) as f:
        records = f.readlines()
    for line in records:
        line = re.sub('\d+:', '', line)
        array = line.strip().split() if line.strip() != '' else None
        encodings.append(array[1:])
        labels.append(int(array[0]))
    return np.array(encodings).astype(float), np.array(labels).astype(int)

def read_tsv(file):
    encodings = []
    labels = []
    with open(file, ) as f:
        records = f.readlines()
    for line in records:
        array = line.strip().split("\t") if line.strip() != '' else None
        encodings.append(array[1:])
        labels.append(int(array[0]))
    return np.array(encodings).astype(float), np.array(labels).astype(int)

def AAC(fastas):
    AA = 'ACDEFGHIKLMNPQRSTVWY'
    encodings = []
    header = ['#']
    for i in AA:
        header.append(i)
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        name = str.split(name, '|')[1]
        count = Counter(sequence)
        for key in count:
            print(count[key])
            count[key] = count[key] / len(sequence)
            print(len(sequence))
        code = [name]
        for aa in AA:
            code.append(count[aa])
        encodings.append(code)
    return encodings

def findMarkedLine(file, win_size):
    # read the original sequence for contact extraction
    originals_sequence = readFasta('/data/residues.fasta')
    originals_sequence = np.array(originals_sequence)

    # mapping ids:
    mapped = pd.read_csv('/data/id_mapping.csv', header=0)
    mapped = np.asarray(mapped)

    # read the no-contact files of certain win_size:
    no_contact = pd.read_table('/data/mapped_segments_length_'+str(win_size)+'.txt', sep="\t", header=None)
    no_contact = np.asarray(no_contact)

    # find the corresponding sequence of content file and original file:
    map_id = str.split(file, '_')[0]
    map_line = np.asarray(mapped[mapped[:, 0] == map_id])
    if len(map_line) == 0:
        print("empty in file " + str(file))
        # continue

    residue_id = map_line[0, 1]
    residue = originals_sequence[originals_sequence[:, 0] == residue_id][:, 1][0]
    # here the marked_line might be a table, as it may contains several sites
    marked_line = np.asarray(no_contact[no_contact[:, 2] == residue_id])
    return residue, marked_line

def savetsv(encodings, file='encoding.tsv'):
    with open(file, 'w') as f:
        if encodings == 0:
            f.write('Descriptor calculation failed.')
        else:
            for i in range(len(encodings[0]) - 1):
                f.write(encodings[0][i] + '\t')
            f.write(encodings[0][-1] + '\n')
            for i in encodings[1:]:
                f.write(i[0] + '\t')
                for j in range(1, len(i) - 1):
                    f.write(str(float(i[j])) + '\t')
                f.write(str(float(i[len(i) - 1])) + '\n')
    return None

def plot_roc_ind(data, out, label_column=0, score_column=2):
    fprIndep, tprIndep, thresholdsIndep = roc_curve(data[:, label_column], data[:, score_column])
    ind_auc = auc(fprIndep, tprIndep)
    fig = plt.figure(0)
    plt.plot(fprIndep, tprIndep, lw=2, alpha=0.7, color='red',
             label='ROC curve (area = %0.2f)' % ind_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc="lower right")
    plt.savefig(out)
    plt.close(0)
    return ind_auc

def calculate_metrics(labels, scores, cutoff=0.5, po_label=1):
    my_metrics = {
        'Sensitivity': 'NA',
        'Specificity': 'NA',
        'Accuracy': 'NA',
        'MCC': 'NA',
        'Recall': 'NA',
        'Precision': 'NA',
        'F1-score': 'NA',
        'Cutoff': cutoff,
    }

    tp, tn, fp, fn = 0, 0, 0, 0
    for i in range(len(scores)):
        if labels[i] == po_label:
            if scores[i] >= cutoff:
                tp = tp + 1
            else:
                fn = fn + 1
        else:
            if scores[i] < cutoff:
                tn = tn + 1
            else:
                fp = fp + 1

    my_metrics['Sensitivity'] = tp / (tp + fn) if (tp + fn) != 0 else 'NA'
    my_metrics['Specificity'] = tn / (fp + tn) if (fp + tn) != 0 else 'NA'
    my_metrics['Accuracy'] = (tp + tn) / (tp + fn + tn + fp)
    my_metrics['MCC'] = (tp * tn - fp * fn) / math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) if (tp + fp) * (
        tp + fn) * (tn + fp) * (tn + fn) != 0 else 'NA'
    my_metrics['Precision'] = tp / (tp + fp) if (tp + fp) != 0 else 'NA'
    my_metrics['Recall'] = my_metrics['Sensitivity']
    my_metrics['F1-score'] = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) != 0 else 'NA'
    return my_metrics

def save_prediction_metrics_ind(m_dict, out):
    with open(out, 'w') as f:
        f.write('#')
        for key in m_dict:
            f.write('\t%s' %key)
        f.write('\n')
        f.write('Indep')
        for key in m_dict:
            f.write('\t%s' %m_dict[key])
        f.write('\n')
    return None