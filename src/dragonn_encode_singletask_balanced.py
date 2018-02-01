# train basset model on single task and multi task predictions
# import functions from my repo
import socket
import sys
hn = socket.gethostname()
if hn == 'aspire':
    dehancerSrc = '/home/ben/projects/dehancer/src/'
elif hn == 'nandi':
    dehancerSrc = '/users/bsiranos/projects/dehancer/src/'
elif hn.startswith('sherlock') or getenv('SHERLOCK') == '1' or getenv('SHERLOCK') == '2':
    dehancerSrc = '/home/bsiranos/projects/dehancer/src/'
else: 
    print('Warning: set location of dehancer repository src file in this script')
    dehancerSrc = '~/projects/dehancer/src/'
sys.path.append(dehancerSrc)
from basset_functions import *
from os.path import join, exists
from os import makedirs, chdir
import numpy as np
from dragonn.models import SequenceDNN
from dragonn.utils import one_hot_encode, get_motif_scores, reverse_complement, encode_fasta_sequences
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18

# info about different tasks
tasks = np.array([i for i in range(32)])
with open('/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/calls_train.txt', 'rU') as inf:
    taskHeader = inf.readline().strip().split('\t')

# model settings
weighted=False
balanced=True

for task in tasks:
    # task = 0
    taskName = taskHeader[task]
    print('Starting task ' + taskName)

    # input files
    if balanced:
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test_balanced/task_' + str(task)
        saveDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/dragonn/' + taskName + '_balanced'
    else:
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
        saveDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/dragonn/' + taskName + '_unbalanced'

    inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa',
                 'calls_train.txt','calls_valid.txt','calls_test.txt']
    inFiles = [join(inFileBase, i) for i in inFileAdd]

    # create save dir if it doesnt exist
    if not exists(saveDir):
        print('making ' + saveDir)
        makedirs(saveDir)

    nb_trainSamples = int(sum(1 for line in open(inFiles[3])) -1)
    nb_validSamples = int(sum(1 for line in open(inFiles[4])) -1)
    nb_testSamples = int(sum(1 for line in open(inFiles[5])) -1)

    print('Training samples: ' + str(nb_trainSamples))
    print('Validation samples: ' + str(nb_validSamples))
    print('Test samples: ' + str(nb_testSamples))
    batchSizeTrain = min(10000, nb_trainSamples)
    batchSizeValid = min(10000, nb_validSamples)
    batchSizeTest = min(10000, nb_testSamples)

    trainGenerator = data_generator(inFiles[0], inFiles[3],
        batchSize=batchSizeTrain, subsetFeatureList=task)
    validGenerator = data_generator(inFiles[1], inFiles[4], 
        batchSize=batchSizeValid, subsetFeatureList=task)
    testGenerator = data_generator(inFiles[2], inFiles[5], 
        batchSize=batchSizeTest, subsetFeatureList=task)


    # make save dir
    if not  exists(saveDir):
        print('making output directory: ' + saveDir)
        makedirs(saveDir) 
    chdir(saveDir)

    # single task to start
    do_hyperparameter_search = False
    seq_length = 2000
    num_hyperparameter_trials = 50
    num_epochs = 100
    use_deep_CNN = False
    use_RNN = False

    trainData = trainGenerator.next()
    validData = validGenerator.next()
    testData = testGenerator.next()
    X_train = trainData[0]
    y_train = trainData[1]
    X_valid = validData[0]
    y_valid = validData[1]
    X_test = testData[0]
    y_test = testData[1]
    # transpose from basset arch to whats expected for dragonn
    X_train = np.transpose(X_train, (0,2,1,3))
    X_valid = np.transpose(X_valid, (0,2,1,3))
    X_test = np.transpose(X_test, (0,2,1,3))
    y_train = np.reshape(y_train, newshape=(len(y_train), 1))
    y_valid = np.reshape(y_valid, newshape=(len(y_valid), 1))
    y_test = np.reshape(y_test, newshape=(len(y_test), 1))

    # Build and train model
    if not do_hyperparameter_search:
        hyperparameters = {'seq_length': seq_length, 'use_RNN': use_RNN,
                           'num_filters': (45,), 'pool_width': 25, 'conv_width': (10,),
                           'L1': 0, 'dropout': 0.2, 'num_epochs': num_epochs}
        if use_deep_CNN:
            hyperparameters.update({'num_filters': (45, 50, 50), 'conv_width': (10, 8, 5)})
        if use_RNN:
            hyperparameters.update({'GRU_size': 35, 'TDD_size': 45})
        model = SequenceDNN(**hyperparameters)
        model.train(X_train, y_train, validation_data=(X_valid, y_valid),
                    save_best_model_to_prefix='best_model')

    else:
        print('Starting hyperparameter search...')
        from dragonn.hyperparameter_search import HyperparameterSearcher, RandomSearch
        fixed_hyperparameters = {'seq_length': seq_length, 'use_RNN': use_RNN, 'num_epochs': num_epochs}
        grid = {'num_filters': ((5, 100),), 'pool_width': (5, 40),
                'conv_width': ((6, 20),), 'dropout': (0, 0.5)}
        if use_deep_CNN:
            grid.update({'num_filters': ((5, 100), (5, 100), (5, 100)),
                         'conv_width': ((6, 20), (6, 20), (6, 20))})
        if use_RNN:
            grid.update({'GRU_size': (10, 50), 'TDD_size': (20, 60)})

        # Backend is RandomSearch; if using Python 2, can also specify MOESearch
        # (requires separate installation)
        searcher = HyperparameterSearcher(SequenceDNN, fixed_hyperparameters, grid, X_train, y_train,
                                          validation_data=(X_valid, y_valid), backend=RandomSearch)
        searcher.search(num_hyperparameter_trials)
        print('Best hyperparameters: {}'.format(searcher.best_hyperparameters))
        model = searcher.best_model


    # Test model
    testResultFile = join(saveDir, 'test_results.txt')
    testResults = model.test(X_test, y_test)
    with open(testResultFile, 'w') as rf:
        rf.write('Test results: {}'.format(testResults) + '\n')
    print('Test results: {}'.format(testResults))


# # Plot DeepLift and ISM scores for the first 10 test examples, and model architecture
# if sys.version[0] == 2:
#     model.plot_deeplift(X_test[:10], output_directory='deeplift_plots')
# model.plot_in_silico_mutagenesis(X_test[:10], output_directory='ISM_plots')
# model.plot_architecture(output_file='architecture_plot.png')

# # save model 
# model.save('test.model')
# model.load('test.model.arch.json')