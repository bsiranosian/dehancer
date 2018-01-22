from os.path import join, exists
from os import makedirs, chdir
import numpy as np
from dragonn.models import SequenceDNN
from dragonn.utils import one_hot_encode, get_motif_scores, reverse_complement, encode_fasta_sequences
try:
    from sklearn.model_selection import train_test_split  # sklearn >= 0.18
except ImportError:
    from sklearn.cross_validation import train_test_split  # sklearn < 0.18
import sys
 
import socket
hn = socket.gethostname()
if hn=='nandi':
    # for use on Nandi 
    fastaDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta'
    saveDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/dragonn/Dnd41_short'
    THEANO_FLAGS='device=cuda3,cuda.root=/usr/bin,floatX=float32,force_device=True'
elif hn=='aspire':
    # for use on aspire
    fastaDir='/home/ben/ak_local/enhancer_conservation/encode_data/broad/fasta'
    saveDir='/home/ben/ak_local/enhancer_conservation/encode_data/broad/dragonn/Dnd41_short'

# make save dir
if not  exists(saveDir):
    print('making output directory: ' + saveDir)
    makedirs(saveDir) 
chdir(saveDir)

# start with one cell line for now 
cell='Dnd41'
do_hyperparameter_search = False
seq_length = 2000
num_hyperparameter_trials = 50
num_epochs = 100
use_deep_CNN = False
use_RNN = False

enhancerFasta=join(fastaDir, cell+'_enhancers_windows_short.fa')
promoterFasta=join(fastaDir, cell+'_promoters_windows_short.fa')
negativeFasta=join(fastaDir, cell+'_negative_windows_short.fa')

# load up the sequences and one hot encode them 
eoh = encode_fasta_sequences(enhancerFasta)
poh = encode_fasta_sequences(promoterFasta)
noh = encode_fasta_sequences(negativeFasta)

# first lets just do enhancer predicition
encodedSequences = np.append(eoh,noh, axis=0)
labels = np.append(np.repeat(True, np.shape(eoh)[0]), np.repeat(False, np.shape(noh)[0]))
labels = np.reshape(labels, newshape=(len(labels), 1))

# split into training and validation
testFraction = 0.2
validationFraction = 0.2
X_train, X_test, y_train, y_test = train_test_split(encodedSequences, labels, test_size=testFraction)
X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, test_size=validationFraction)

# add reverse complement, to training set only
X_train = np.concatenate((X_train, reverse_complement(X_train)))
y_train = np.concatenate((y_train, y_train))

# Randomly splitting data into training and test 
random_order = np.arange(len(X_train))
np.random.shuffle(random_order)
X_train = X_train[random_order]
y_train = y_train[random_order]

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
print('Test results: {}'.format(model.test(X_test, y_test)))
# Plot DeepLift and ISM scores for the first 10 test examples, and model architecture
if sys.version[0] == 2:
    model.plot_deeplift(X_test[:10], output_directory='deeplift_plots')
model.plot_in_silico_mutagenesis(X_test[:10], output_directory='ISM_plots')
model.plot_architecture(output_file='architecture_plot.png')

# save model 
model.save('test.model')
model.load('test.model.arch.json')