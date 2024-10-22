import pandas as pd
import glob
import os
import os.path as op
from collections import Counter
import sys
import numpy as np
from sklearn import neighbors
import pandas as pd
from math import log

from .tools import safe_makedir



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


def merge_all(phage_df_loc, recruit_df_loc, outfile, training_file='/mnt/scgc_nfs/opt/viruscope/virus-training.csv'):
    ''' merge recruitment and phage orf dataframes, calculate probabilities and write to outfile
    Args:
    '''
    phage_df = pd.read_csv(phage_df_loc)
    recruit_df = pd.read_csv(recruit_df_loc)
    merged = recruit_df.merge(phage_df.reset_index(), on='contig', how='outer')
    result = merged.merge(virus_class(merged, training_file=training_file), on='contig', how='outer')
    result.to_csv(outfile, index=False)
    return result


def write_batch_summaries(contig_dir, wd, training_file='/mnt/scgc_nfs/opt/viruscope/virus-training.csv'):
    ''' summarize results and calculate statistics at the end of the viruSCope run.
    Args:
        contig_dir(str): path to directory containing contigs in fasta format
        wd (str): path to viruscope main output directory
    Returns:
        path (str): path to directory containing viruscope summary tables
    '''

    summary_dir = safe_makedir(op.join(wd, "summary"))
    contigs = glob.glob(op.join(contig_dir, "*.f*"))
    _blast = lambda wd, name: op.join(wd, 'blast', '{}_blast_summary.csv'.format(name))
    _diamond = lambda wd, name: op.join(wd, 'diamond', '{}_proteins_mg_diamond_recruitment_tbl.csv'.format(name))
    _summary = lambda wd, name: op.join(wd, 'summary', '{}_summary.csv'.format(name))

    for c in contigs:
        name = op.basename(c).split(".")[0].split("_")[0]
        if op.exists(_summary(wd, name)):
            print(_summary(wd, name), "exists!  Skipping.", file=sys.stdout)
            continue

        try:
            out = merge_all(_blast(wd, name), _diamond(wd, name), _summary(wd, name), training_file=training_file)
        except Exception as inst:
            print('WARNING, unable to create summary table for {name}: {inst}'.format(name=name, inst=inst), file=sys.stdout)
            continue

    print("viruscope summary tables written to {}".format(summary_dir), file=sys.stdout)
    return summary_dir
