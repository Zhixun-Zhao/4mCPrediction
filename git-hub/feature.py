import numpy as np
from repDNA.psenac import PseDNC


# sequence to 4-bit binary
def binary(seq):
    vector = ''
    for i in seq:
        if i == 'A':
            vector += '1 0 0 0 '
        elif i == 'C':
            vector += '0 1 0 0 '
        elif i == 'G':
            vector += '0 0 1 0 '
        else:
            vector += '0 0 0 1 '
    vector = vector.strip()
    return vector


# overlapping chemical and physical features
# def opf(seq):
#     vector = ''
#     count = 0
#     count_a = 0
#     count_c = 0
#     count_g = 0
#     count_t = 0
#     for i in seq:
#         count += 1
#         if i == 'A':
#             count_a += 1
#             density = count_a / count
#             vector += '1 1 1' + ' ' + str(density) + ' '
#         elif i == 'C':
#             count_c += 1
#             density = count_c / count
#             vector += '0 0 1' + ' ' + str(density) + ' '
#         elif i == 'G':
#             count_g += 1
#             density = count_g / count
#             vector += '1 0 0' + ' ' + str(density) + ' '
#         else:
#             count_t += 1
#             density = count_t / count
#             vector += '0 1 0' + ' ' + str(density) + ' '
#     vector = vector.strip()
#     return vector

def density(seq):
    vector = ''
    count = 0
    count_a = 0
    count_c = 0
    count_g = 0
    count_t = 0
    for i in seq:
        count += 1
        if i == 'A':
            count_a += 1
            density = count_a / count
            vector += str(density) + ' '
        elif i == 'C':
            count_c += 1
            density = count_c / count
            vector += str(density) + ' '
        elif i == 'G':
            count_g += 1
            density = count_g / count
            vector += str(density) + ' '
        else:
            count_t += 1
            density = count_t / count
            vector += str(density) + ' '
    vector = vector.strip()
    return vector


# k-mer
def kmer(seq):
    count_1mer = np.zeros((1, 4))
    for i in seq:
        if i == 'A':
            count_1mer[0,0] += 1
        elif i == 'C':
            count_1mer[0,1] += 1
        elif i == 'G':
            count_1mer[0,2] += 1
        else:
            count_1mer[0,3] += 1
    density_1mer = count_1mer / len(seq)

    count_2mer = np.zeros((1, 16))
    for i in range(len(seq) - 1):
        if seq[i:i + 2] == 'AA':
            count_2mer[0,0] += 1
        if seq[i:i + 2] == 'AC':
            count_2mer[0,1] += 1
        if seq[i:i + 2] == 'AG':
            count_2mer[0,2] += 1
        if seq[i:i + 2] == 'AT':
            count_2mer[0,3] += 1
        if seq[i:i + 2] == 'CA':
            count_2mer[0,4] += 1
        if seq[i:i + 2] == 'CC':
            count_2mer[0,5] += 1
        if seq[i:i + 2] == 'CG':
            count_2mer[0,6] += 1
        if seq[i:i + 2] == 'CT':
            count_2mer[0,7] += 1
        if seq[i:i + 2] == 'GA':
            count_2mer[0,8] += 1
        if seq[i:i + 2] == 'GC':
            count_2mer[0,9] += 1
        if seq[i:i + 2] == 'GG':
            count_2mer[0,10] += 1
        if seq[i:i + 2] == 'GT':
            count_2mer[0,11] += 1
        if seq[i:i + 2] == 'TA':
            count_2mer[0,12] += 1
        if seq[i:i + 2] == 'TC':
            count_2mer[0,13] += 1
        if seq[i:i + 2] == 'TG':
            count_2mer[0,14] += 1
        if seq[i:i + 2] == 'TT':
            count_2mer[0,15] += 1
    density_2mer = count_2mer / (len(seq) - 1)
    vector = np.hstack((density_1mer, density_2mer))
    vector = str(vector).replace('\n', '').replace('[', '').replace(']', '')
    return vector


# KSNPFS
def ksnpfs(seq, k):
    count_2mer = np.zeros((1, 16))
    for i in range(len(seq) - k - 1):
        if seq[i] == 'A' and seq[i + k + 1] == 'A':
            count_2mer[0, 0] += 1
        if seq[i] == 'A' and seq[i + k + 1] == 'C':
            count_2mer[0, 1] += 1
        if seq[i] == 'A' and seq[i + k + 1] == 'G':
            count_2mer[0, 2] += 1
        if seq[i] == 'A' and seq[i + k + 1] == 'T':
            count_2mer[0, 3] += 1
        if seq[i] == 'C' and seq[i + k + 1] == 'A':
            count_2mer[0, 4] += 1
        if seq[i] == 'C' and seq[i + k + 1] == 'C':
            count_2mer[0, 5] += 1
        if seq[i] == 'C' and seq[i + k + 1] == 'G':
            count_2mer[0, 6] += 1
        if seq[i] == 'C' and seq[i + k + 1] == 'T':
            count_2mer[0, 7] += 1
        if seq[i] == 'G' and seq[i + k + 1] == 'A':
            count_2mer[0, 8] += 1
        if seq[i] == 'G' and seq[i + k + 1] == 'C':
            count_2mer[0, 9] += 1
        if seq[i] == 'G' and seq[i + k + 1] == 'G':
            count_2mer[0, 10] += 1
        if seq[i] == 'G' and seq[i + k + 1] == 'T':
            count_2mer[0, 11] += 1
        if seq[i] == 'T' and seq[i + k + 1] == 'A':
            count_2mer[0, 12] += 1
        if seq[i] == 'T' and seq[i + k + 1] == 'C':
            count_2mer[0, 13] += 1
        if seq[i] == 'T' and seq[i + k + 1] == 'G':
            count_2mer[0, 14] += 1
        if seq[i] == 'T' and seq[i + k + 1] == 'T':
            count_2mer[0, 15] += 1
    density_2mer = count_2mer / (len(seq) - k - 1)
    density_2mer = str(density_2mer).replace('\n', '').replace('[', '').replace(']', '')

    return density_2mer


#PseDNC
def psednc(seq):
    psednc = PseDNC(lamada=3, w=0.05)
    vector = psednc.make_psednc_vec([seq])
    return str(vector).replace('[', '').replace(']', '').replace(',', '')


#Hidden Markov Chain
def hmc_line_prob(seq):
    count = np.zeros((len(seq)-1, 16))
    for i in range(len(seq)-1):
        if seq[i] == 'A' and seq[i + 1] == 'A':
            count[i, 0] += 1
        if seq[i] == 'A' and seq[i + 1] == 'C':
            count[i, 1] += 1
        if seq[i] == 'A' and seq[i + 1] == 'G':
            count[i, 2] += 1
        if seq[i] == 'A' and seq[i + 1] == 'T':
            count[i, 3] += 1
        if seq[i] == 'C' and seq[i + 1] == 'A':
            count[i, 4] += 1
        if seq[i] == 'C' and seq[i  + 1] == 'C':
            count[i, 5] += 1
        if seq[i] == 'C' and seq[i + 1] == 'G':
            count[i, 6] += 1
        if seq[i] == 'C' and seq[i + 1] == 'T':
            count[i, 7] += 1
        if seq[i] == 'G' and seq[i + 1] == 'A':
            count[i, 8] += 1
        if seq[i] == 'G' and seq[i + 1] == 'C':
            count[i, 9] += 1
        if seq[i] == 'G' and seq[i + 1] == 'G':
            count[i, 10] += 1
        if seq[i] == 'G' and seq[i + 1] == 'T':
            count[i, 11] += 1
        if seq[i] == 'T' and seq[i + 1] == 'A':
            count[i, 12] += 1
        if seq[i] == 'T' and seq[i + 1] == 'C':
            count[i, 13] += 1
        if seq[i] == 'T' and seq[i + 1] == 'G':
            count[i, 14] += 1
        if seq[i] == 'T' and seq[i + 1] == 'T':
            count[i, 15] += 1
    return count


def hmm(train, test):
    count_all = np.zeros((len(train[0].strip())-1, 16))
    for line in train:
        line = line.strip()
        count = hmc_line_prob(line)
        count_all = np.add(count_all, count)

    for k in range(len(count_all)):
        count_1 = count_all[k, 1] + count_all[k, 2] + count_all[k, 3] + count_all[k, 0]
        count_2 = count_all[k, 5] + count_all[k, 6] + count_all[k, 7] + count_all[k, 4]
        count_3 = count_all[k, 9] + count_all[k, 10] + count_all[k, 11] + count_all[k, 8]
        count_4 = count_all[k, 13] + count_all[k, 14] + count_all[k, 15] + count_all[k, 12]
        for t in range(4):
            count_all[k, t] = count_all[k, t] / count_1
        for t in range(4, 8):
            count_all[k, t] = count_all[k, t] / count_2
        for t in range(8, 12):
            count_all[k, t] = count_all[k, t] / count_3
        for t in range(12, 16):
            count_all[k, t] = count_all[k, t] / count_4

    test_sample = test[0].strip()
    test_data = np.zeros((len(test), len(test_sample)-1))
    for j in range(len(test)):
        seq = test[j].strip()
        for i in range(len(seq)-1):
            if seq[i] == 'A' and seq[i + 1] == 'A':
                test_data[j, i] = count_all[i, 0]
            if seq[i] == 'A' and seq[i + 1] == 'C':
                test_data[j, i] = count_all[i, 1]
            if seq[i] == 'A' and seq[i + 1] == 'G':
                test_data[j, i] = count_all[i, 2]
            if seq[i] == 'A' and seq[i + 1] == 'T':
                test_data[j, i] = count_all[i, 3]
            if seq[i] == 'C' and seq[i + 1] == 'A':
                test_data[j, i] = count_all[i, 4]
            if seq[i] == 'C' and seq[i  + 1] == 'C':
                test_data[j, i] = count_all[i, 5]
            if seq[i] == 'C' and seq[i + 1] == 'G':
                test_data[j, i] = count_all[i, 6]
            if seq[i] == 'C' and seq[i + 1] == 'T':
                test_data[j, i] = count_all[i, 7]
            if seq[i] == 'G' and seq[i + 1] == 'A':
                test_data[j, i] = count_all[i, 8]
            if seq[i] == 'G' and seq[i + 1] == 'C':
                test_data[j, i] = count_all[i, 9]
            if seq[i] == 'G' and seq[i + 1] == 'G':
                test_data[j, i] = count_all[i, 10]
            if seq[i] == 'G' and seq[i + 1] == 'T':
                test_data[j, i] = count_all[i, 11]
            if seq[i] == 'T' and seq[i + 1] == 'A':
                test_data[j, i] = count_all[i, 12]
            if seq[i] == 'T' and seq[i + 1] == 'C':
                test_data[j, i] = count_all[i, 13]
            if seq[i] == 'T' and seq[i + 1] == 'G':
                test_data[j, i] = count_all[i, 14]
            if seq[i] == 'T' and seq[i + 1] == 'T':
                test_data[j, i] = count_all[i, 15]

    return test_data

