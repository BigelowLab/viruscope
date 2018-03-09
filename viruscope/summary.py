import pandas as pd
import glob
import os
import os.path as op
from collections import Counter

from bisect import bisect_left
import numpy as np
from sklearn import neighbors
import pandas as pd
from math import log




def virus_class(vsdf, test_cols = ['viral_phage_gene_fraction','ratio_virus_bacteria'],
                training_file = '/mnt/scgc_nfs/opt/viruscope/virus-training.csv',
                k = 27, score_col='score'):
    ''' function for determining virus class using KNeighborsClassifier
    Args:
        vsdf (pandas.DataFrame): viruscope dataframe having the the two calculated test columns and a contigs column included
        test_cols (list): list of names for test columns to be used for KNeighbors Classifier
        training_file (path): path to the training file for the knn algorithm
        k (int): number of nearest neighbors
        score_col(str): name of score column to use for training
    Returns:
        dataframe with columns: contig, virus_class, probability
    '''
    def _pval_by_count(clf, test_set, training_df,
                  test_cols=['viral_phage_gene_fraction','ratio_virus_bacteria'],
                 score_col='score'):

        '''
        function for calculating probability that matches the previously used R function

        Args:
            clf: KNeighborsClassifier object
            test_set: dataframe used to test knn algorithm on
            training_df: dataframe used to train knn
            test_cols: columns in above dataframes used to train the knn algorithm
            score_col: column in training algorithm that reports the score
        Returns:
            list of pvalues in same order as rows in test_set dataframe
        '''

        ns = clf.kneighbors(test[test_cols], return_distance=False)
        pvals = []

        for l in ns:
            ndf = training_df.loc[l,:]
            pvals.append(Counter(ndf['score']).most_common()[0][1] / len(ndf))
        return pvals

    train_in = pd.read_table(training_file, sep = ',').dropna(axis = 0, how = 'any').reset_index()
    train_in.rename(columns={'ratio.virus.bacteria':'ratio_virus_bacteria'}, inplace=True)
    train_in[test_cols[0]] = np.log10(train_in[test_cols[0]] + 0.001)
    train_in[test_cols[1]] = np.log10(train_in[test_cols[1]] + 0.001)
    train_label = list(train_in.iloc[:,-1].values)
    in_train = train_in[test_cols]
    in_train.columns = test_cols

    # make the classifier and fit the training the data

    clf = neighbors.KNeighborsClassifier(k, weights='distance', algorithm='auto')
    clf.fit(in_train, train_label)

    test = vsdf[['contig'] + test_cols].dropna()
    if len(test) == 0:
        return pd.DataFrame(columns=['contig','virus_class','virus_prob','virus_prob_new'])
    test[test_cols[0]] = np.log10(test[test_cols[0]] + 0.001)
    test[test_cols[1]] = np.log10(test[test_cols[1]] + 0.001)
    try:
        test['virus_class'] = clf.predict(test[test_cols])
    except Exception as inst:
        print(test)
        raise inst
    probs = clf.predict_proba(test[test_cols])
    # built in probability function
    test['prob 0'] = [i[0] for i in probs]
    test['prob 1'] = [i[1] for i in probs]
    test['virus_prob_new'] = [l['prob 0'] if l['virus_class'] == 0 else l['prob 1'] for i, l in test.iterrows()]

    test['virus_prob'] = _pval_by_count(clf, test, train_in, test_cols, score_col)

    return test[['contig','virus_class','virus_prob','virus_prob_new']]


def merge_all(phage_df, recruit_df, outfile):
    merged = recruit_df.merge(phage_df.reset_index(), on='contig', how='outer')
    result = merged.merge(virus_class(merged), on='contig', how='outer')
    result.to_csv(outfile, index=False)
    return result
    
