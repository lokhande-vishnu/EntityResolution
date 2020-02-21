#!/usr/bin/python
# -*- coding: utf-8 -*-

from future.builtins import next

import os
import csv
import re
import logging
import optparse
import collections
import numpy as np
from scipy import stats

import dedupe
from unidecode import unidecode

from utils import *

import scipy.io

optp = optparse.OptionParser()
optp.add_option('-v', '--verbose', dest='verbose', action='count',
                help='Increase verbosity (specify multiple times for more)'
                )
(opts, args) = optp.parse_args()
log_level = logging.WARNING 
if opts.verbose:
    if opts.verbose == 1:
        log_level = logging.INFO
    elif opts.verbose >= 2:
        log_level = logging.DEBUG
logging.getLogger().setLevel(log_level)

# ## Setup
exp = 'results_1/'
seed = 103
USE_CHECKPOINT = False
input_file = 'csv_example_messy_input.csv'
output_file = exp + 'VAL_csv_example_output.csv'
zeus_output_file = exp + 'VAL_ZEUS_csv_example_output.csv'
settings_file = exp + 'csv_example_learned_settings'
golden_file = 'csv_example_input_with_true_ids.csv'
val_golden_file = exp + 'VAL_csv_example_golden.csv'
F_file = exp + 'F_sample'

if not os.path.exists(exp):
    os.mkdir(exp)
    print("Directory " , exp ,  " Created ")
else:
    print("Directory " , exp ,  " already exists")

'''
Transforms pairwise observation scores into a graph data structure.
'''
def makeGraph(scores):
    F = {'pairwise':[], 'neigh':[], 'mapO2N':{}, 'mapN2O':{}}

    # Create Index Maps
    row_id = 0
    for score in scores:
        n1 = score[0][0]
        n2 = score[0][1]
        if n1 not in F['mapO2N']:
            F['mapO2N'][n1] = row_id + 1
            F['mapN2O'][row_id + 1] = n1
            row_id += 1
        if n2 not in F['mapO2N']:
            F['mapO2N'][n2] = row_id + 1
            F['mapN2O'][row_id + 1] = n2
            row_id += 1
    
    # Pairwise Costs
    min_n = float('inf')
    for score in scores:
        n1 = score[0][0]
        n2 = score[0][1]

        # Transform probability to a score between -1 to 1.
        w = -2*(score[1] - 0.5)
        if F['mapO2N'][n1] < F['mapO2N'][n2]:
            F['pairwise'].append([F['mapO2N'][n1], F['mapO2N'][n2], w])
        else:
            F['pairwise'].append([F['mapO2N'][n2], F['mapO2N'][n1], w])

        if min_n > min(F['mapO2N'][n1], F['mapO2N'][n2]):
            min_n = min(F['mapO2N'][n1], F['mapO2N'][n2])

    if min_n != 1:
        print("ERROR: index should start with 1")
        exit(1)

    # Neigh

    ## Strategy 1
    for score in scores:
        n1 = score[0][0]
        n2 = score[0][1]
        if F['mapO2N'][n1] < F['mapO2N'][n2]:
            F['neigh'].append([F['mapO2N'][n1], F['mapO2N'][n2]])
        else:
            F['neigh'].append([F['mapO2N'][n2], F['mapO2N'][n1]])

    scipy.io.savemat(F_file + '.mat', F)
    np.save(F_file + '.npy' , F)

'''
Print the statistics of the true cluster sizes within the dataset.
''' 
def printGraphStats(data_d, label_d):
    dupe_d = collections.defaultdict(list)

    for cid in label_d.values():
        if cid not in dupe_d:
            dupe_d[cid] = 1.
        else:
            dupe_d[cid] += 1.

    if 'x' in dupe_d :
        del dupe_d['x']

    dupe_d_values = list(dupe_d.values())
    print('Graph Stats ...')
    print('Number of clusters:', len(dupe_d_values))
    print('Mean size of cluster:', np.mean(dupe_d_values))
    print('Median size of cluster:', np.median(dupe_d_values))
    print('Mode size of cluster:', stats.mode(dupe_d_values))
    print('MAX size of cluster:', max(dupe_d_values))
    print('MIN size of cluster:', min(dupe_d_values))

'''
Save observations along with a cluster to which it belongs to.
'''        
def saveCluster(singleton_id, data_d, cluster_membership, canonical_keys, output_file, input_file):
    with open(output_file, 'w') as f_output, open(input_file) as f_input:
        writer = csv.writer(f_output)
        reader = csv.reader(f_input)

        heading_row = next(reader)
        heading_row.insert(0, 'confidence_score')
        heading_row.insert(0, 'Cluster ID')
        for key in canonical_keys:
            heading_row.append('canonical_' + key)

        writer.writerow(heading_row)

        for row in reader:
            row_id = int(row[0])
            if row_id in cluster_membership:
                cluster_id = cluster_membership[row_id]["cluster id"]
                canonical_rep = cluster_membership[row_id]["canonical representation"]
                row.insert(0, cluster_membership[row_id]['confidence'])
                row.insert(0, cluster_id)
                for key in canonical_keys:
                    row.append(canonical_rep[key].encode('utf8'))
                writer.writerow(row)
            elif row_id in data_d:
                row.insert(0, None)
                row.insert(0, singleton_id)
                singleton_id += 1
                for key in canonical_keys:
                    row.append(None)
                writer.writerow(row)

'''
Splits a given dataset into training set and validation set
'''
def split(data_d, label_d, split_ratio = 0.5, seed = 42):
    train_data_d, train_label_d = {}, {}
    val_data_d, val_label_d = {}, {}
    np.random.seed(seed)
    coins = np.random.binomial(1, split_ratio, len(data_d))
    i = -1
    for rid in data_d:
        i += 1
        if coins[i] > 0:
            train_data_d[rid] = data_d[rid]
            train_label_d[rid] = label_d[rid]
        else:
            val_data_d[rid] = data_d[rid]
            val_label_d[rid] = label_d[rid]
    return train_data_d, train_label_d, val_data_d, val_label_d

'''
Samples a subset from a given dataset for training
'''
def sampleTrainData(data_d, label_d, train_ratio = 0.5, seed = 42):
    train_data_d, train_label_d = {}, {}
    val_data_d, val_label_d = {}, {}
    np.random.seed(seed)
    coins = np.random.binomial(1, train_ratio, len(data_d))
    i = -1
    for rid in data_d:
        i += 1
        if coins[i] > 0:
            train_data_d[rid] = data_d[rid]
            train_label_d[rid] = label_d[rid]
            
        val_data_d[rid] = data_d[rid]
        val_label_d[rid] = label_d[rid]
    return train_data_d, train_label_d, val_data_d, val_label_d

'''
Trains a classifer on a small subset using the true labels
'''                
def train(data_d, label_d, sample_ratio = 1.0):
    # Define the fields dedupe will pay attention to
    fields = [
        {'field' : 'Site name', 'type': 'String'},
        {'field' : 'Address', 'type': 'String'},
        {'field' : 'Zip', 'type': 'Exact', 'has missing' : True},
        {'field' : 'Phone', 'type': 'String', 'has missing' : True},
        ]

    print('start labeling...')
    training_pairs = makeTrainingPairs(data_d, label_d, sample_ratio, seed)

    print('match: ', len(training_pairs['match']),
          'distinct:', len(training_pairs['distinct']))

    # Create a new deduper object and pass our data model to it.
    deduper = dedupe.Dedupe(fields)

    ## To train dedupe, we feed it a sample of records.
    deduper.sample(data_d, len(data_d)*len(data_d))

    
    print('marking pairs ...')
    deduper.markPairs(training_pairs)
                            
    print('training ...')
    deduper.train()

    with open(settings_file, 'wb') as sf:
        deduper.writeSettings(sf)


    return deduper

'''
Assess the performance on the validation dataset
'''        
def test(data_d):
    # BLOCK
    print('blocking ...')
    blocked_pairs = deduper.blockData(data_d)
    print(' ')

    # SCORE
    print('scoring ...')
    scores = deduper.scoreBlocks(blocked_pairs)

    return scores


print('importing data ...')
data_d, label_d, gold_fieldnames = readGoldData(golden_file)

print('graph stats ...')
print(' ')
printGraphStats(data_d, label_d)

print('split data ...')
train_data_d, train_label_d, val_data_d, val_label_d = sampleTrainData(data_d, label_d,
                                                                       train_ratio = 0.5,
                                                                       seed = seed)

print(len(train_data_d), len(train_label_d))
print(len(val_data_d), len(val_label_d))

with open(val_golden_file, 'w') as g_output:
    writer = csv.DictWriter(g_output, fieldnames=gold_fieldnames)
    writer.writeheader()
    for rid, row in val_data_d.items():
        row['True Id'] = label_d[rid]
        writer.writerow(row)
                            

# If a settings file already exists, we'll just load that and skip training
print('begin training ...')
if USE_CHECKPOINT and os.path.exists(settings_file):
    print('reading from', settings_file)
    with open(settings_file, 'rb') as f:
        deduper = dedupe.StaticDedupe(f)
else:
    deduper = train(train_data_d, train_label_d, 0.01)


val_scores = test(val_data_d)

#######################################################################
 
# Make Graph for ZeuS      
makeGraph(val_scores)

# HIERARCHICAL CLUSTERING
print('hierarchical clustering ...')
val_clustered_dupes = list(deduper.clusterScores(val_scores)) #, threshold))
cluster_membership = {}
cluster_id = 0
for (cluster_id, cluster) in enumerate(val_clustered_dupes):
    id_set, scores = cluster
    cluster_d = [val_data_d[c] for c in id_set]
    canonical_rep = dedupe.canonicalize(cluster_d)
    for record_id, score in zip(id_set, scores):
        cluster_membership[record_id] = {
            "cluster id" : cluster_id,
            "canonical representation" : canonical_rep,
            "confidence": score
        }
singleton_id = cluster_id + 1
saveCluster(singleton_id,
            val_data_d,
            cluster_membership,
            canonical_rep.keys(),
            output_file, input_file)



