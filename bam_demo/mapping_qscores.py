import numpy as np
import pandas as pd

# bam_demo otuput
df_map = pd.read_csv("Ch7_scRNA_AVITI_qscore.txt", header=None)
qscore_map = dict()
for _, row in df_map.iterrows():
    if row[2] not in qscore_map:
        qscore_map[int(row[2])] = [int(row[0])]
    else:
        qscore_map[int(row[2])].append(int(row[0]))

# get_qc output
df_orig = pd.read_csv("quality_score_csvs/Ch7_scRNA_AVITI_R2_QC.csv", header=None)
count_orig = df_orig.to_numpy()

count_new = np.zeros(df_orig.shape, dtype=np.int64)
for qscore in qscore_map:
    count_new[:, qscore] = count_orig[:, qscore_map[qscore]].sum(axis=1)

df_new = pd.DataFrame(count_new)
df_new.to_csv("ch7_aviti_qc.corrected.txt", header=False, index=False)
