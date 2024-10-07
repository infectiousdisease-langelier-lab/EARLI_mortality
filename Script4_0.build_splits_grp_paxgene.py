'''

'''

# load libraries and functions
# 2.2.2
import pandas as pd
# 1.4.2
from sklearn.model_selection import StratifiedKFold

# set paths
data_path = "Inputs/"
results_path = "Outputs/"

# set params
target_cat = 'Hospital_Death'
groups_to_include = ['G1', 'G2']
suffix_ = ''.join(groups_to_include)

suffix_ = suffix_ + '_' + target_cat

results_path = results_path

group_dict = {'G1': '1_Sepsis+BldCx+',
              'G2': '2_Sepsis+OtherCx+',
              'G4': '4_NO_Sepsis'}

#########################
# DATA LOADING
#########################
# counts data
paxgene_cnt_data = pd.read_csv(data_path + "earli_counts_kallisto_mortality.csv", index_col=0)

# meta data
paxgene_meta_data = pd.read_csv(data_path + "Cleaned_Metadata_070324.csv")
paxgene_meta_data = paxgene_meta_data.sort_values(by=['Barcode'])

#########################
# DATA PRE-PROCESSING
#########################
# samples of interest (sepsis) only
paxgene_meta_data['EARLI_Barcode'] = ['EARLI_' + str(x) for x in paxgene_meta_data.Barcode]
paxgene_meta_data = paxgene_meta_data.loc[paxgene_meta_data.Group.isin([group_dict[x] for x in groups_to_include]), :]

paxgene_samples = list(set(paxgene_meta_data.EARLI_Barcode.values) & set(paxgene_cnt_data.columns.values))
paxgene_samples = sorted(paxgene_samples)

paxgene_cnt_data = paxgene_cnt_data[paxgene_samples]
paxgene_cnt_data = paxgene_cnt_data.T

paxgene_meta_data.index = paxgene_meta_data.EARLI_Barcode
paxgene_meta_data = paxgene_meta_data.loc[paxgene_samples, :]

pd.DataFrame(data=[paxgene_cnt_data.index, paxgene_meta_data.Group, paxgene_meta_data[target_cat]]).T.to_csv(results_path + "included_paxgene_kallisto_samples_" + suffix_ + ".csv")

#########################
# CREATE FOLDS
#########################
# split data to have a held-out set
skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=123)

i = 0
for train_index, test_index in skf.split(paxgene_cnt_data, paxgene_meta_data[target_cat].values):
        X_train, X_test = paxgene_cnt_data.iloc[train_index, :], paxgene_cnt_data.iloc[test_index, :]
        y_train, y_test = paxgene_meta_data[target_cat].values[train_index], paxgene_meta_data[target_cat].values[test_index]

        pd.DataFrame(data=X_test.index).to_csv(results_path + "test_labels_" + str(i + 1) + "_" + suffix_ + ".csv")
        i = i+1

