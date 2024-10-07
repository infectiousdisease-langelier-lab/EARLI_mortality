'''

'''

# load libraries and functions
# to load custom functions
from custom_classes import CustomBaggingClassifier
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

# 1.4.2
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC

# set paths
data_path = "Inputs/"
results_path = "Outputs/"

# set params
# "create" from scratch or "load" existing dump files
mode_ = "create"

# outcome
target_cat = 'Hospital_Death'
groups_to_include = ['G1', 'G2']
suffix_ = ''.join(groups_to_include)

suffix_ = suffix_ + "_" + target_cat

group_dict = {'G1': '1_Sepsis+BldCx+',
              'G2': '2_Sepsis+OtherCx+',
              'G3': '3_Sepsis+Cx-',
              'G4': '4_NO_Sepsis',
              'G5': '5_Unclear'}

excl_filt_ntrpm = False
incl_filt_ntrpm = False
incl_exl_filt_ntrpm = True

#########################
# DATA LOADING
#########################
# counts data
if excl_filt_ntrpm:
    bgc_data = pd.read_csv(results_path + "ntrpm_bgc_data_cast_exclude_list.csv", index_col=0)
else:
    if incl_filt_ntrpm:
        bgc_data = pd.read_csv(results_path + "ntrpm_bgc_data_cast_include_list.csv", index_col=0)
    else:
        if incl_exl_filt_ntrpm:
            bgc_data = pd.read_csv(results_path + "ntrpm_bgc_data_cast_include_exclude_list.csv", index_col=0)
            bgc_data_G4 = pd.read_csv(results_path + "ntrpm_bgc_data_cast_include_exclude_list_G4.csv", index_col=0)
        else:
            bgc_data = pd.read_csv(data_path + "bgc_data_cast.csv", index_col=0)

# meta data
paxgene_meta_data = pd.read_csv(data_path + "Cleaned_Metadata_070324.csv")

#########################
# DATA PRE-PROCESSING
#########################
bgc_data = bgc_data.T
bgc_data_G4 = bgc_data_G4.T

bgc_data_G4 = np.log(bgc_data_G4+1)

# samples of interest (sepsis) only
paxgene_meta_data_G4 = paxgene_meta_data.loc[paxgene_meta_data.Group.isin([group_dict[x] for x in ['G4']]), :]
paxgene_meta_data = paxgene_meta_data.loc[paxgene_meta_data.Group.isin([group_dict[x] for x in groups_to_include]), :]

paxgene_meta_data.index = ['EARLI_' + str(x) for x in paxgene_meta_data.Barcode]
paxgene_meta_data_G4.index = ['EARLI_' + str(x) for x in paxgene_meta_data_G4.Barcode]

paxgene_samples = list(set(paxgene_meta_data.index.values) & set(bgc_data.index.values))
paxgene_samples = sorted(paxgene_samples)

bgc_data = bgc_data.loc[paxgene_samples, :]
paxgene_meta_data = paxgene_meta_data.loc[paxgene_samples, :]

vars_to_adjust = []
cov_to_add = []

C_val = 1 # 1 for full frame
n_est = 100 # 100 for full frame
max_feat = 0.9 # 0.5 for full frame
max_samp = 1.0 # 1.0 for full frame

n_genes = 50

model_ = "bsvm"

feat_sel = 'dgea' #dgea lasso

suffix_dgea = suffix_

suffix_out = feat_sel + "_" + str(C_val) + "_" + str(n_est) + "_" + str(max_feat) + "_" + str(max_samp) + "_" + str(n_genes)

if excl_filt_ntrpm:
    suffix_out = suffix_out + "_excl_filt_ntrpm"
if incl_filt_ntrpm:
    suffix_out = suffix_out + "_incl_filt_ntrpm"
if incl_exl_filt_ntrpm:
    suffix_out = suffix_out + "_incl_excl_filt_ntrpm"

original_set = True
if original_set:
    suffix_out = suffix_out + "_original_set"

#########################
# FIT MODEL
#########################
res_dict_bsvm = []
cvs_aucroc_test_bsvm = []
cvs_aucroc_test_bsvm_G4 = []

for i in range(0, 10):
    i = i+1
    i = str(i)

    test_ids = pd.read_csv(results_path + "test_labels_" + i + "_" + suffix_ + ".csv", index_col=0)
    
    #########################
    # DATA PRE-PROCESSING
    #########################
    # prepare data
    # split data to have a held-out set
    bgc_data_train = bgc_data.loc[bgc_data.index.isin(test_ids.iloc[:, 0].values)==False, :]

    bgc_data_train = bgc_data_train.loc[:, bgc_data_train.astype(bool).sum(axis=0)>2]

    bgc_data_test = bgc_data.loc[bgc_data.index.isin(test_ids.iloc[:, 0].values)==True, bgc_data_train.columns.values]

    bgc_data_train = np.log(bgc_data_train+1)
    bgc_data_test = np.log(bgc_data_test+1)

    paxgene_target_cat_train = paxgene_meta_data.loc[bgc_data_train.index.values, target_cat]
    paxgene_target_cat_test = paxgene_meta_data.loc[bgc_data_test.index.values, target_cat]

    #########################
    # BagSVM - CV within a CV for grid search and RFE as part of a pipeline
    # n_jobs should be set to 1/None (or at least not -1 in both CustomRFECV and GridSearchCV)
    # to avoid nested parallelism
    #########################
    # evaluate performance using roc_auc values
    # transform, scale, filter, classifier
    pipe_bsvm = Pipeline([ ('scale', StandardScaler()),
                                    ('bsvmc', CustomBaggingClassifier(estimator=LinearSVC(max_iter=100000, 
                                                                                               C=C_val, dual=True),
                                                                    random_state=123,
                                                                    max_features=max_feat,
                                                                    max_samples=max_samp,
                                                                    n_estimators=n_est
                                                                    ))])

    search_bsvm = pipe_bsvm

    # fit model
    if mode_ == "create":
        search_bsvm.fit(bgc_data_train, paxgene_target_cat_train)
        dump(search_bsvm, results_path + "bgc_" + model_ + "_dump_" + suffix_out + "_" + i + "_" + suffix_dgea + ".joblib")
    else: # + "" +
        if mode_ == "load":
            search_bsvm = load(results_path + "bgc_" + model_ + "_dump_" + suffix_out + "_" + i + "_" + suffix_dgea + ".joblib")

    # test on full test set
    probs_bsvm = search_bsvm.predict_proba(bgc_data_test)[:, 1]
    roc_auc_bsvm = roc_auc_score(paxgene_target_cat_test, probs_bsvm)

    print(bgc_data_train.shape)
    print(roc_auc_score(paxgene_target_cat_train, search_bsvm.predict_proba(bgc_data_train)[:, 1]))
    print(roc_auc_bsvm)

    # feature importance
    best_vars = np.concatenate([bgc_data_train.columns.values[x] for x in search_bsvm.named_steps['bsvmc'].estimators_features_]).ravel()
    best_coefs = np.concatenate([search_bsvm.named_steps['bsvmc'].estimators_[x].coef_ for x in range(0, len(search_bsvm.named_steps['bsvmc'].estimators_))]).ravel()

    pd.DataFrame({'best_vars':best_vars, 'best_coefs':best_coefs}).to_csv(results_path + "top_vars_" + model_ + "_" + suffix_out + "_" + i +"_" + suffix_dgea + ".csv")

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
    
    # G4
    bgc_data_G4_fold = bgc_data_G4.loc[:, bgc_data_train.columns.values]
    probs_bsvm_G4 = search_bsvm.predict_proba(bgc_data_G4_fold)[:, 1]
    roc_auc_bsvm_G4 = roc_auc_score(paxgene_meta_data_G4.loc[bgc_data_G4_fold.index.values, target_cat], probs_bsvm_G4)

    cvs_aucroc_test_bsvm_G4.append({"labels": paxgene_meta_data_G4.loc[bgc_data_G4_fold.index.values, target_cat],
                                 "probs": probs_bsvm_G4,
                                 "classes": search_bsvm.classes_[1],
                                 "roc_auc": roc_auc_bsvm_G4})

# save to CSV
pd.DataFrame(data=res_dict_bsvm).to_csv(results_path + "bgc_summary_table_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

cvs_aucroc_test_all_bsvm = {"samples":[y for z in [x["labels"].index.values for x in cvs_aucroc_test_bsvm] for y in z],
                            "y_true":[y for z in [x["labels"] for x in cvs_aucroc_test_bsvm] for y in z],
                            "y_pred_prob":[y for z in[x["probs"] for x in cvs_aucroc_test_bsvm] for y in z],
                            "class_pos":[y for z in[np.repeat(x["classes"], len(x["labels"])) for x in cvs_aucroc_test_bsvm] for y in z]}

pd.DataFrame.from_dict(cvs_aucroc_test_all_bsvm).to_csv(results_path + "bgc_pred_table_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

# G4
cvs_aucroc_test_all_bsvm_G4 = {"samples":[y for z in [x["labels"].index.values for x in cvs_aucroc_test_bsvm_G4] for y in z],
                            "y_true":[y for z in [x["labels"] for x in cvs_aucroc_test_bsvm_G4] for y in z],
                            "y_pred_prob":[y for z in[x["probs"] for x in cvs_aucroc_test_bsvm_G4] for y in z],
                            "class_pos":[y for z in[np.repeat(x["classes"], len(x["labels"])) for x in cvs_aucroc_test_bsvm_G4] for y in z]}

pd.DataFrame.from_dict(cvs_aucroc_test_all_bsvm_G4).to_csv(results_path + "bgc_pred_table_G4_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".csv")

print(str(np.mean([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm_G4]).round(2)))
print(str(np.std([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm_G4]).round(2)))

#########################
# AUC-ROC curves
#########################
# bsvm
# plot ROC curve for the testing set
tprs_ = []
base_fpr = np.linspace(0, 1, 101)

plt.figure()
for i in range(0, 10):
    cv_id = str(i + 1)
    fpr_test, tpr_test, _ = roc_curve(cvs_aucroc_test_bsvm[i]["labels"],
                                      cvs_aucroc_test_bsvm[i]["probs"],
                                      pos_label=cvs_aucroc_test_bsvm[i]["classes"])

    tpr_ = interp(base_fpr, fpr_test, tpr_test)
    tpr_[0] = 0.0
    tprs_.append(tpr_)

tprs_ = np.array(tprs_)
mean_tprs = tprs_.mean(axis=0)
std_tprs = tprs_.std(axis=0)

tprs_upper = np.minimum(mean_tprs + std_tprs, 1)
tprs_lower = mean_tprs - std_tprs

plt.plot(base_fpr, mean_tprs, 'red', label='Cross-validation splits, AUC=' +
                                           str(np.mean([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm]).round(2)) + " (" +
                                           str(np.std([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm]).round(2)) + ")")
plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='red', alpha=0.2)

plt.plot([0, 1], [0, 1], color='grey', linewidth=1)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic (ROC) Curve - Test sets')
plt.legend()

plt.savefig(results_path + "bgc_test_auc_roc_" + model_ + "_" + suffix_out + "_" + suffix_dgea + ".pdf")

print(np.mean([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm]).round(2))
print(np.std([x["roc_auc"].round(2) for x in cvs_aucroc_test_bsvm]).round(2))
