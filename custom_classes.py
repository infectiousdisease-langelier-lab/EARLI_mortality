'''
custom_classes.py
created on Nov 2 2020
lucile.neyton@ucsf.edu

This file contains different classes necessary to run the script in build_classifiers.py.
More specifically, it contains a transformer class allowing to transform RNA-Seq counts data using the variance
stabilising transformation.
A modified version of the BaggingClassifier (sklearn, 0.23.2) class is provided as well to permit the extraction
of features necessary to rank variables given their importance.
'''

# load libraries and functions
# 1.4.2
from joblib import Parallel, delayed, effective_n_jobs
# 1.26.4
import numpy as np
# 3.5.16
from rpy2 import robjects
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri
# 1.4.2
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import BaggingClassifier

# activate R objects and load libraries
numpy2ri.activate()
pandas2ri.activate()
importr("DESeq2")


# define classes and functions
class VstTransformer(BaseEstimator, TransformerMixin):
    def __init__(self):
        self.disp_train = None

    def fit(self, X, *_):

        # compute the distributions to be used for future transformations
        robjects.globalenv["cnt_data_train"] = X.T
        robjects.globalenv["meta_data_train"] = robjects.r('data.frame(dummy_column = '
                                                           'c(rep(1,dim(cnt_data_train)[2])), row.names=colnames(cnt_data_train))')
        
        robjects.globalenv["dds_train"] = robjects.r('DESeqDataSetFromMatrix(countData = cnt_data_train,'
                                                     'colData = meta_data_train, design = ~1)')
        robjects.globalenv["dds_train"] = robjects.r('estimateSizeFactors(dds_train)')
        robjects.globalenv["dds_train"] = robjects.r('estimateDispersions(dds_train)')
        self.disp_train = robjects.r('dispersionFunction(dds_train)')

        return self

    def transform(self, X, *_):
        print(X.shape)
        robjects.globalenv["cnt_data"] = X.T
        robjects.globalenv["meta_data"] = robjects.r('data.frame(dummy_column = c(rep(1,dim(cnt_data)[2])), row.names=colnames(cnt_data))')
        robjects.globalenv["dds"] = robjects.r('DESeqDataSetFromMatrix(countData = cnt_data,'
                                               'colData = meta_data, design = ~1)')
        robjects.globalenv["dds"] = robjects.r('estimateSizeFactors(dds)')
        robjects.globalenv["dispersionFunction(dds)"] = self.disp_train
        X = robjects.r('t(assay(varianceStabilizingTransformation(dds, blind = FALSE)))')
        return X



# define classes and functions
class VstTransformerWithCovariates(BaseEstimator, TransformerMixin):
    def __init__(self, cov_name):
        self.disp_train = None
        self.cov_name = cov_name

    def fit(self, X, *_):
        # compute the distributions to be used for future transformations
        robjects.globalenv["cnt_data_train"] = X.loc[:, X.columns.values != self.cov_name].T
        robjects.globalenv["meta_data_train"] = robjects.r('data.frame(dummy_column = '
                                                           'c(rep(1,dim(cnt_data_train)[2])), row.names=colnames(cnt_data_train))')
        robjects.globalenv["dds_train"] = robjects.r('DESeqDataSetFromMatrix(countData = cnt_data_train,'
                                                     'colData = meta_data_train, design = ~1)')
        robjects.globalenv["dds_train"] = robjects.r('estimateSizeFactors(dds_train)')
        robjects.globalenv["dds_train"] = robjects.r('estimateDispersions(dds_train)')
        self.disp_train = robjects.r('dispersionFunction(dds_train)')

        return self

    def transform(self, X, *_):
        print(X.shape)
        robjects.globalenv["covs"] = X.loc[:, self.cov_name]

        cov_save = X.loc[:, self.cov_name]

        robjects.globalenv["cnt_data"] = X.loc[:, X.columns.values != self.cov_name].T
        robjects.globalenv["meta_data"] = robjects.r('data.frame(dummy_column = c(rep(1,dim(cnt_data)[2])), row.names=colnames(cnt_data))')
        robjects.globalenv["dds"] = robjects.r('DESeqDataSetFromMatrix(countData = cnt_data,'
                                               'colData = meta_data, design = ~1)')
        robjects.globalenv["dds"] = robjects.r('estimateSizeFactors(dds)')
        robjects.globalenv["dispersionFunction(dds)"] = self.disp_train
        X = robjects.r('t(assay(varianceStabilizingTransformation(dds, blind = FALSE)))')

        X = np.append(X, cov_save.to_numpy().reshape((X.shape[0], 1)), axis=1)

        return X

from sklearn.utils import check_random_state, indices_to_mask
from sklearn.utils.validation import _check_sample_weight, has_fit_parameter
import numbers
from sklearn.utils.random import sample_without_replacement
from warnings import warn
import itertools
MAX_INT = np.iinfo(np.int32).max

def _partition_estimators(n_estimators, n_jobs):
    """Private function used to partition estimators between jobs."""
    # Compute the number of jobs
    n_jobs = min(effective_n_jobs(n_jobs), n_estimators)

    # Partition estimators between jobs
    n_estimators_per_job = np.full(n_jobs, n_estimators // n_jobs,
                                   dtype=int)
    n_estimators_per_job[:n_estimators % n_jobs] += 1
    starts = np.cumsum(n_estimators_per_job)

    return n_jobs, n_estimators_per_job.tolist(), [0] + starts.tolist()

def _parallel_build_estimators(n_estimators, ensemble, X, y, sample_weight,
                               seeds, total_n_estimators, verbose):
    """Private function used to build a batch of estimators within a job."""
    # Retrieve settings
    n_samples, n_features = X.shape
    max_features = ensemble._max_features
    max_samples = ensemble._max_samples
    bootstrap = ensemble.bootstrap
    bootstrap_features = ensemble.bootstrap_features
    support_sample_weight = has_fit_parameter(ensemble.estimator_,
                                              "sample_weight")
    if not support_sample_weight and sample_weight is not None:
        raise ValueError("The base estimator doesn't support sample weight")

    # Build estimators
    estimators = []
    estimators_features = []

    for i in range(n_estimators):
        if verbose > 1:
            print("Building estimator %d of %d for this parallel run "
                  "(total %d)..." % (i + 1, n_estimators, total_n_estimators))

        random_state = seeds[i]
        estimator = ensemble._make_estimator(append=False,
                                             random_state=random_state)

        # Draw random feature, sample indices
        features, indices = _generate_bagging_indices(random_state,
                                                      bootstrap_features,
                                                      bootstrap, n_features,
                                                      n_samples, max_features,
                                                      max_samples)

        ############################ EDIT START ###########################################
        # all samples drawn from the same class
        while len(list(set(y[indices])))==1:
            print(y[indices])
            seeds[i] = seeds[i] + 1
            random_state = random_state + 1
            estimator = ensemble._make_estimator(append=False,
                                                 random_state=random_state)
            features, indices = _generate_bagging_indices(random_state,
                                                          bootstrap_features,
                                                          bootstrap, n_features,
                                                          n_samples, max_features,
                                                          max_samples)
        ############################ EDIT END #############################################

        # Draw samples, using sample weights, and then fit
        if support_sample_weight:
            if sample_weight is None:
                curr_sample_weight = np.ones((n_samples,))
            else:
                curr_sample_weight = sample_weight.copy()

            if bootstrap:
                sample_counts = np.bincount(indices, minlength=n_samples)
                curr_sample_weight *= sample_counts
            else:
                not_indices_mask = ~indices_to_mask(indices, n_samples)
                curr_sample_weight[not_indices_mask] = 0

            estimator.fit(X[:, features], y, sample_weight=curr_sample_weight)

        else:
            estimator.fit((X[indices])[:, features], y[indices])

        estimators.append(estimator)
        estimators_features.append(features)

    return estimators, estimators_features, seeds


def _generate_bagging_indices(random_state, bootstrap_features,
                              bootstrap_samples, n_features, n_samples,
                              max_features, max_samples):
    """Randomly draw feature and sample indices."""
    # Get valid random state
    random_state = check_random_state(random_state)

    # Draw indices
    feature_indices = _generate_indices(random_state, bootstrap_features,
                                        n_features, max_features)
    sample_indices = _generate_indices(random_state, bootstrap_samples,
                                       n_samples, max_samples)

    return feature_indices, sample_indices


def _generate_indices(random_state, bootstrap, n_population, n_samples):
    """Draw randomly sampled indices."""
    # Draw sample indices
    if bootstrap:
        indices = random_state.randint(0, n_population, n_samples)
    else:
        indices = sample_without_replacement(n_population, n_samples,
                                             random_state=random_state)

    return indices


# @property decorator used so that these two properties can be extracted without having to call a getter
class CustomBaggingClassifier(BaggingClassifier):
    @property
    def feature_importances_(self):
        return np.mean(np.power([mod.coef_ for mod in self.estimators_], 2), axis=0)

    def _fit(self, X, y, max_samples=None, max_depth=None, sample_weight=None):
        random_state = check_random_state(self.random_state)

        # Convert data (X is required to be 2d and indexable)
        X, y = self._validate_data(
            X, y, accept_sparse=['csr', 'csc'], dtype=None,
            force_all_finite=False, multi_output=True
        )
        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X, dtype=None)

        # Remap output
        n_samples, self.n_features_ = X.shape
        self._n_samples = n_samples
        y = self._validate_y(y)

        # Check parameters
        self._validate_estimator()

        if max_depth is not None:
            self.estimator_.max_depth = max_depth

        # Validate max_samples
        if max_samples is None:
            max_samples = self.max_samples
        elif not isinstance(max_samples, numbers.Integral):
            max_samples = int(max_samples * X.shape[0])

        if not (0 < max_samples <= X.shape[0]):
            raise ValueError("max_samples must be in (0, n_samples]")

        # Store validated integer row sampling value
        self._max_samples = max_samples

        # Validate max_features
        if isinstance(self.max_features, numbers.Integral):
            max_features = self.max_features
        elif isinstance(self.max_features, float):
            max_features = self.max_features * self.n_features_
        else:
            raise ValueError("max_features must be int or float")

        if not (0 < max_features <= self.n_features_):
            raise ValueError("max_features must be in (0, n_features]")

        max_features = max(1, int(max_features))

        # Store validated integer feature sampling value
        self._max_features = max_features

        # Other checks
        if not self.bootstrap and self.oob_score:
            raise ValueError("Out of bag estimation only available"
                             " if bootstrap=True")

        if self.warm_start and self.oob_score:
            raise ValueError("Out of bag estimate only available"
                             " if warm_start=False")

        if hasattr(self, "oob_score_") and self.warm_start:
            del self.oob_score_

        if not self.warm_start or not hasattr(self, 'estimators_'):
            # Free allocated memory, if any
            self.estimators_ = []
            self.estimators_features_ = []

        n_more_estimators = self.n_estimators - len(self.estimators_)

        if n_more_estimators < 0:
            raise ValueError('n_estimators=%d must be larger or equal to '
                             'len(estimators_)=%d when warm_start==True'
                             % (self.n_estimators, len(self.estimators_)))

        elif n_more_estimators == 0:
            warn("Warm-start fitting without increasing n_estimators does not "
                 "fit new trees.")
            return self

        # Parallel loop
        n_jobs, n_estimators, starts = _partition_estimators(n_more_estimators,
                                                             self.n_jobs)
        total_n_estimators = sum(n_estimators)

        # Advance random state to state after training
        # the first n_estimators
        if self.warm_start and len(self.estimators_) > 0:
            random_state.randint(MAX_INT, size=len(self.estimators_))

        seeds = random_state.randint(MAX_INT, size=n_more_estimators)
        self._seeds = seeds

        all_results = Parallel(n_jobs=n_jobs, verbose=self.verbose,
                               **self._parallel_args())(
            delayed(_parallel_build_estimators)(
                n_estimators[i],
                self,
                X,
                y,
                sample_weight,
                seeds[starts[i]:starts[i + 1]],
                total_n_estimators,
                verbose=self.verbose)
            for i in range(n_jobs))

        # Reduce
        self.estimators_ += list(itertools.chain.from_iterable(
            t[0] for t in all_results))
        self.estimators_features_ += list(itertools.chain.from_iterable(
            t[1] for t in all_results))
        self._seeds += list(itertools.chain.from_iterable(
            t[2] for t in all_results))

        if self.oob_score:
            self._set_oob_score(X, y)

        return self


class DGEA_filter(BaseEstimator, TransformerMixin):
    def __init__(self, vars_to_keep):
        self.vars_to_keep = vars_to_keep

    def fit(self, X, y=None):
        return self

    def transform(self, X, y=None):
        print(X[:, self.vars_to_keep].shape)
        return X[:, self.vars_to_keep]

