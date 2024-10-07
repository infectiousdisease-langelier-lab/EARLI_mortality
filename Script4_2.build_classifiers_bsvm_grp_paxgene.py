
# load libraries and functions
from custom_classes import VstTransformer, VstTransformerWithCovariates, CustomBaggingClassifier, DGEA_filter
# 1.4.2
from joblib import dump, load
# 3.9.1
import matplotlib.pyplot as plt
# 1.26.4
import numpy as np
from numpy import interp
# 2.2.2
import pandas as pd
#
import re
# 3.5.16
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri
numpy2ri.activate()
pandas2ri.activate()
importr("DESeq2")
# 0.13.2
import seaborn as sns
# 1.4.2
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

# set paths
data_path = "Inputs/"
results_path = "Outputs/"

# set params
# "create" from scratch or "load" existing dump files
mode_ = "load"

# outcome
target_cat = 'Hospital_Death'
groups_to_include = ['G1', 'G2']
suffix_ = ''.join(groups_to_include)

suffix_ = suffix_ + "_" + target_cat

group_dict = {'G1': '1_Sepsis+BldCx+',
              'G2': '2_Sepsis+OtherCx+',
              'G4': '4_NO_Sepsis'}

vars_to_adjust = []
#vars_to_adjust = ["apacheiii_scaled", "microbial_mass"]

add_cov = False
cov_to_add = 'staph_or_hhv5' #'microbial_mass'

if add_cov:
    vars_to_adjust = []
else:
    cov_to_add = []

C_val = 0.1
n_est = 100
n_genes = 50
max_feat = 0.1 #0.25
max_samp = 1.0 #1.0

model_ = "bsvm"

feat_sel = 'dgea'

suffix_dgea = suffix_
for var_to_adjust in vars_to_adjust:
    suffix_dgea = suffix_dgea + "_" + var_to_adjust + "_AND"

suffix_dgea = re.sub('_AND$', '', suffix_dgea)

suffix_out = feat_sel + "_" + str(C_val) + "_" + str(n_est) + "_" + str(max_feat) + "_" + str(max_samp) + "_" + str(n_genes)
if add_cov:
    suffix_out = suffix_out + "_addedCov_" + cov_to_add

# counts data
paxgene_full_cnt_data = pd.read_csv(results_path + "paxgene_cnt_data.csv", index_col=0)
paxgene_full_cnt_data = paxgene_full_cnt_data.T

# paxgene meta data (to order main metadata frame)
paxgene_full_meta_data = pd.read_csv(results_path + "paxgene_meta_data.csv", index_col=0)

# G4
paxgene_g4_meta_data = paxgene_full_meta_data.loc[paxgene_full_meta_data.Group.isin([group_dict[x] for x in ['G4']]), :]
paxgene_g4_cnt_data = paxgene_full_cnt_data.loc[paxgene_g4_meta_data.index, :]

# filter for groups of interest
paxgene_full_meta_data = paxgene_full_meta_data.loc[paxgene_full_meta_data.Group.isin([group_dict[x] for x in groups_to_include]), :]
paxgene_full_cnt_data = paxgene_full_cnt_data.loc[paxgene_full_meta_data.index, :]

#########################
# FIT MODEL
#########################
res_dict_bsvm = []
cvs_aucroc_test_bsvm = []
res_dict_bsvm_g4 = []
cvs_aucroc_test_bsvm_g4 = []

for i in range(0, 10):
    i = i+1
    i = str(i)

    #########################
    # DATA LOADING
    #########################
    # counts data
    paxgene_cnt_data = pd.read_csv(data_path + "paxgene_cnt_data_" + i + "_" + suffix_ + ".csv", index_col=0)
    paxgene_cnt_data = paxgene_cnt_data.T

    # paxgene meta data (to order main metadata frame)
    paxgene_meta_data = pd.read_csv(data_path + "paxgene_meta_data_" + i + "_" + suffix_ + ".csv", index_col=0)

    # DGEA list
    sig_res_paxgene = pd.read_csv(results_path + "sig_res_paxgene_" + i + "_" + suffix_dgea + ".csv", index_col=0)

    sig_res_paxgene = sig_res_paxgene.iloc[0:n_genes, :]

    paxgene_full_cnt_data = paxgene_full_cnt_data.drop("EARLI_11942", errors="ignore")
    paxgene_cnt_data = paxgene_cnt_data.drop("EARLI_11942", errors="ignore")
    paxgene_meta_data = paxgene_meta_data.drop("EARLI_11942", errors="ignore")
    paxgene_full_meta_data = paxgene_full_meta_data.drop("EARLI_11942", errors="ignore")

    # prepare data
    # split data to have a held-out set
    paxgene_cnt_data_train = paxgene_cnt_data
    paxgene_cnt_data_test = paxgene_full_cnt_data.loc[
                            paxgene_full_cnt_data.index.isin(paxgene_cnt_data_train.index.values) == False,
                            paxgene_full_cnt_data.columns.isin(paxgene_cnt_data_train.columns.values)]

    paxgene_target_cat_train = paxgene_meta_data.loc[:, target_cat]
    paxgene_target_cat_test = paxgene_full_meta_data.loc[paxgene_cnt_data_test.index, target_cat]

    #########################
    # DATA PRE-PROCESSING
    #########################
    # list genes used as input
    name_vars = paxgene_cnt_data.columns.values

    if feat_sel == 'dgea':
        y = sig_res_paxgene.index.values
        

    y = [name_vars.tolist().index(x) for x in y]

    #########################
    # BagSVM - CV within a CV for grid search and RFE as part of a pipeline
    # n_jobs should be set to 1/None (or at least not -1 in both CustomRFECV and GridSearchCV)
    # to avoid nested parallelism
    #########################
    if add_cov:
        paxgene_meta_data['staph_or_hhv5'] = (paxgene_meta_data.staph | paxgene_meta_data.hhv5)
        paxgene_full_meta_data['staph_or_hhv5'] = (paxgene_full_meta_data.staph | paxgene_full_meta_data.hhv5)

        paxgene_cnt_data_train[cov_to_add] = paxgene_meta_data.loc[:, cov_to_add]
        paxgene_cnt_data_test[cov_to_add] = paxgene_full_meta_data.loc[paxgene_cnt_data_test.index, cov_to_add]
        y = np.append(y, len(name_vars))

    # evaluate performance using roc_auc values
    # transform, scale, filter, classifier
    if add_cov:
        if model_ == "bsvm":
            pipe_bsvm = Pipeline([('norm', VstTransformerWithCovariates(cov_name=cov_to_add)),
                                    ('scale', StandardScaler()),
                                    ('filt', DGEA_filter(vars_to_keep=y)),
                                    ('bsvmc', CustomBaggingClassifier(estimator=LinearSVC(max_iter=10000, C=C_val),
                                                                    random_state=123,
                                                                    max_features=max_feat,
                                                                    max_samples=max_samp,
                                                                    n_estimators=n_est
                                                                    ))])

    else:
        if model_ == "bsvm":
            pipe_bsvm = Pipeline([('norm', VstTransformer()),
                                    ('scale', StandardScaler()),
                                    ('filt', DGEA_filter(vars_to_keep=y)),
                                    ('bsvmc', CustomBaggingClassifier(estimator=LinearSVC(max_iter=10000, C=C_val),
                                                                    random_state=123,
                                                                    max_features=max_feat,
                                                                    max_samples=max_samp,
                                                                    n_estimators=n_est
                                                                    ))])
                                                                    

    search_bsvm = pipe_bsvm

    # fit model
    if mode_ == "create":
        search_bsvm.fit(paxgene_cnt_data_train, paxgene_target_cat_train)
        dump(search_bsvm, results_path + "paxgene_" + model_ + "_dump_" + suffix_out + "_" + i + "_" + suffix_dgea + ".joblib")
    else: # + "" +
        if mode_ == "load":
            search_bsvm = load(results_path + "paxgene_" + model_ + "_dump_" + suffix_out + "_" + i + "_" + suffix_dgea + ".joblib")

    # test on full test set
    probs_bsvm = search_bsvm.predict_proba(paxgene_cnt_data_test)[:, 1]
    roc_auc_bsvm = roc_auc_score(paxgene_target_cat_test, probs_bsvm)

    print(roc_auc_score(paxgene_target_cat_train, search_bsvm.predict_proba(paxgene_cnt_data_train)[:, 1]))
    print(roc_auc_bsvm)

    #########################
    # Performance summary
    #########################
    # performance summary
    res_dict_bsvm = res_dict_bsvm + [roc_auc_bsvm.round(2)]

    # evaluate on test data
    cvs_aucroc_test_bsvm.append({"labels": paxgene_target_cat_test,
                                    "probs": probs_bsvm,
                                    "classes": search_bsvm.classes_[1],
                                    "roc_auc": roc_auc_bsvm})

    # test on G4
    # test on full test set
    if add_cov:
        paxgene_g4_meta_data['staph_or_hhv5'] = (paxgene_g4_meta_data.staph | paxgene_g4_meta_data.hhv5)
        paxgene_g4_cnt_data[cov_to_add] = paxgene_g4_meta_data.loc[:, cov_to_add]

    probs_g4_bsvm = search_bsvm.predict_proba(paxgene_g4_cnt_data.loc[:, paxgene_g4_cnt_data.columns.isin(paxgene_cnt_data_train.columns.values)])[:, 1]
    roc_auc_g4_bsvm = roc_auc_score(paxgene_g4_meta_data.loc[:, target_cat], probs_g4_bsvm)

    pd.DataFrame.from_dict({"samples": paxgene_g4_meta_data.index.values, "probs": probs_g4_bsvm}).to_csv(
        results_path + "pred_g4_" + i + "_" + suffix_dgea + ".csv")

    # performance summary
    res_dict_bsvm_g4 = res_dict_bsvm_g4 + [roc_auc_g4_bsvm.round(2)]

    # evaluate on test data
    cvs_aucroc_test_bsvm_g4.append({"labels": paxgene_g4_meta_data.loc[:, target_cat],
                                 "probs": probs_g4_bsvm,
                                 "classes": search_bsvm.classes_[1],
                                 "roc_auc": roc_auc_g4_bsvm})

# save to CSV
pd.DataFrame(data=res_dict_bsvm).to_csv(results_path + "summary_table_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

cvs_aucroc_test_all_bsvm = {"samples":[y for z in [x["labels"].index.values for x in cvs_aucroc_test_bsvm] for y in z],
                            "y_true":[y for z in [x["labels"] for x in cvs_aucroc_test_bsvm] for y in z],
                            "y_pred_prob":[y for z in[x["probs"] for x in cvs_aucroc_test_bsvm] for y in z],
                            "class_pos":[y for z in[np.repeat(x["classes"], len(x["labels"])) for x in cvs_aucroc_test_bsvm] for y in z]}

pd.DataFrame.from_dict(cvs_aucroc_test_all_bsvm).to_csv(results_path + "pred_table_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

# G4
pd.DataFrame(data=res_dict_bsvm_g4).to_csv(results_path + "summary_table_g4_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

cvs_aucroc_test_all_bsvm_g4 = {"samples":[y for z in [x["labels"].index.values for x in cvs_aucroc_test_bsvm_g4] for y in z],
                            "y_true":[y for z in [x["labels"] for x in cvs_aucroc_test_bsvm_g4] for y in z],
                            "y_pred_prob":[y for z in[x["probs"] for x in cvs_aucroc_test_bsvm_g4] for y in z],
                            "class_pos":[y for z in[np.repeat(x["classes"], len(x["labels"])) for x in cvs_aucroc_test_bsvm_g4] for y in z]}

pd.DataFrame.from_dict(cvs_aucroc_test_all_bsvm_g4).to_csv(results_path + "pred_table_g4_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")


