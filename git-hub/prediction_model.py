import feature
import joblib

fr = open('.../Test_data.txt', 'r')
fw = open('.../Test_data_feature.txt', 'w')
lines = fr.readlines()
for line in lines:
    if line[0] != '>':
        sequence = line.strip()
        fea1 = feature.binary(sequence)  # 164 bit
        fea2 = feature.density(sequence)  # 40 bit
        fea3 = feature.kmer(sequence)  # 4 + 16 = 20 bit
        fea4 = feature.ksnpfs(sequence, 1)  # 16 bit
        fea5 = feature.ksnpfs(sequence, 2)  # 16 bit
        fea6 = feature.ksnpfs(sequence, 3)  # 16 bit
        fea7 = feature.psednc(sequence)  # 16 + 3 bit
        fw.write(fea1 + " " + fea2 + " " + fea3 + " " + fea4 + " " + fea5 + " " + fea6 + " " + fea7 + "\n")
fw.close()

test_data = np.loadtxt('.../Test_data_feature.txt', dtype=float)
model = joblib.load('./species_name.pkl')
predicted_label = model.predict(test_data_k)

probability_threshold = k
predicted_label[predicted_label > k] = 1
predicted_label[~(predicted_label > k)] = 0