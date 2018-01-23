# using the basset architecture to predict enhancer and promoter sites
# basset arch from anna
import numpy as np 
import pdb
from Bio import SeqIO
from dragonn.utils import one_hot_encode

def data_generator(windowFastaFile, windowCallFile, batchSize=256, mode='train'):
    numGenerated=0
    # get number of fasta sequences in the file
    totalEntries = sum(1 for line in open(windowFastaFile)) /2

    # can get both train and validation data from this
    # set the split fraction here 
    trainFraction = 0.8
    validationFraction = 0.2
    nb_trainSamples = int(np.ceil(batchSize * trainFraction))
    nb_validationSamples = int(batchSize - nb_trainSamples)

    if mode=='train':

        with open(windowFastaFile, 'rU') as ff, open(windowCallFile, 'rU') as cf:
            fastaGen = SeqIO.parse(ff, "fasta")
            # read first header line of calls
            trash=cf.readline()
            # callHeader = np.fromstring(cf.readline(), sep='\t', dtype=str)

            while (numGenerated + batchSize) < totalEntries:
                # keep a set for training
                batchSequences = np.array([str(fastaGen.next().seq.upper()) for i in range(nb_trainSamples)])
                batchOneHot = one_hot_encode(batchSequences)
                batchOneHot = np.transpose(batchOneHot, (0,2,1,3))
                batchCalls = np.array([np.fromstring(cf.readline(), sep='\t',dtype=bool) for i in range(nb_trainSamples)])
                # print(batchCalls.shape)
                # batchCalls = np.reshape(batchCalls, newshape=(batchCalls.shape[0], batchCalls.shape[1], 1))
            
                # trash a set for validation
                trash = [fastaGen.next() for i in range(nb_validationSamples)]
                trash = [cf.readline() for i in range(nb_validationSamples)]
                numGenerated += batchSize
                yield tuple([batchOneHot, batchCalls])

    elif mode=='validation':
        with open(windowFastaFile, 'rU') as ff, open(windowCallFile, 'rU') as cf:

            fastaGen = SeqIO.parse(ff, "fasta")
            # read first header line of calls
            trash=cf.readline()
            # callHeader = np.fromstring(cf.readline(), sep='\t', dtype=str)

            while (numGenerated + batchSize) < totalEntries:
                # trash a set for training
                trash = [fastaGen.next() for i in range(nb_trainSamples)]
                trash = [cf.readline() for i in range(nb_trainSamples)]
                # keep a set for validation
                batchSequences = np.array([str(fastaGen.next().seq.upper()) for i in range(nb_validationSamples)])
                batchOneHot = one_hot_encode(batchSequences)
                batchOneHot = np.transpose(batchOneHot, (0,2,1,3))
                batchCalls = np.array([np.fromstring(cf.readline(), sep='\t',dtype=bool) for i in range(nb_validationSamples)])
                # batchCalls = np.reshape(batchCalls, newshape=(batchCalls.shape[0], batchCalls.shape[1], 1))
            
                numGenerated += batchSize
                yield tuple([batchOneHot, batchCalls])


def getModelGivenModelOptionsAndWeightInits(w0_file,w1_file,init_weights):
    import numpy as np
    np.random.seed(1234)
    import keras;
    from keras.models import Sequential
    from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
    from keras.layers.convolutional import Convolution2D, MaxPooling2D
    from keras.optimizers import Adadelta, SGD, RMSprop;
    # from keras.losses import *;
    from keras.constraints import maxnorm;
    from keras.layers.normalization import BatchNormalization
    from keras.regularizers import l1, l2    
    from keras import backend as K

    model=Sequential()
    model.add(Convolution2D(nb_filter=300,nb_row=1, nb_col=19 ,input_shape=(4,1,2000)))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1,4)))

    model.add(Convolution2D(nb_filter=200,nb_row=1, nb_col=11 ))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1,4)))

    model.add(Convolution2D(nb_filter=200,nb_row=1, nb_col=7 ))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1,4)))

    model.add(Flatten())
    model.add(Dense(1000))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(Dropout(0.3))

    model.add(Dense(1000))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(Dropout(0.3))

    model.add(Dense(32))
    model.add(Activation("sigmoid"))

    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    if (w0_file!=None and w1_file!=None):
        w0=[float(i) for i in open(w0_file,'r').read().strip().split('\n')]
        w1=[float(i) for i in open(w1_file,'r').read().strip().split('\n')]
        loss=get_weighted_binary_crossentropy(w0_weights=w0,w1_weights=w1)
        model.compile(optimizer=adam,loss=loss)
    else:
        model.compile(optimizer=adam,loss="binary_crossentropy")
        # pdb.set_trace()
    return model


# and do same prediction as before
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
    saveDir='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/dragonn/Dnd41_basset_short'
elif hn=='aspire':
    # for use on aspire
    fastaDir='/home/ben/ak_local/enhancer_conservation/encode_data/broad/fasta'
    saveDir='/home/ben/ak_local/enhancer_conservation/encode_data/broad/dragonn/Dnd41_basset_short'

# make save dir
if not  exists(saveDir):
    print('making output directory: ' + saveDir)
    makedirs(saveDir) 
chdir(saveDir)

# # start with one cell line for now 
# cell='Dnd41'
# do_hyperparameter_search = False
# seq_length = 2000
# num_hyperparameter_trials = 50
# num_epochs = 100
# use_deep_CNN = False
# use_RNN = False

# enhancerFasta=join(fastaDir, cell+'_enhancers_windows_short.fa')
# promoterFasta=join(fastaDir, cell+'_promoters_windows_short.fa')
# negativeFasta=join(fastaDir, cell+'_negative_windows_short.fa')

# # load up the sequences and one hot encode them 
# eoh = encode_fasta_sequences(enhancerFasta)
# eoh = np.transpose(eoh, (0,2,1,3))

# poh = encode_fasta_sequences(promoterFasta)
# poh = np.transpose(poh, (0,2,1,3))
# noh = encode_fasta_sequences(negativeFasta)
# noh = np.transpose(noh, (0,2,1,3))

# # first lets just do enhancer predicition
# encodedSequences = np.append(eoh,noh, axis=0)
# labels = np.append(np.repeat(True, np.shape(eoh)[0]), np.repeat(False, np.shape(noh)[0]))
# labels = np.reshape(labels, newshape=(len(labels), 1))

# # split into training and validation
# testFraction = 0.2
# validationFraction = 0.2
# X_train, X_test, y_train, y_test = train_test_split(encodedSequences, labels, test_size=testFraction)
# X_train, X_valid, y_train, y_valid = train_test_split(X_train, y_train, test_size=validationFraction)

# # add reverse complement, to training set only
# X_train = np.concatenate((X_train, reverse_complement(X_train)))
# y_train = np.concatenate((y_train, y_train))

# # Randomly splitting data into training and test 
# random_order = np.arange(len(X_train))
# np.random.shuffle(random_order)
# X_train = X_train[random_order]
# y_train = y_train[random_order]

model=getModelGivenModelOptionsAndWeightInits(None, None, None)
# model.fit(X_train, y_train, validation_data=(X_valid, y_valid))


windowFastaFile = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
windowCallFile = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls/all.calls'
totalEntries = sum(1 for line in open(windowFastaFile)) /2
totalEntries = 10000

trainFraction = 0.8
validationFraction = 0.2
nb_trainSamples = int(np.ceil(totalEntries * trainFraction))
nb_validationSamples = int(totalEntries - nb_trainSamples)

model.fit_generator(generator=data_generator(windowFastaFile, windowCallFile),
                    samples_per_epoch =  nb_trainSamples, nb_epoch=100, show_accuracy=True,
                    validation_data = data_generator(windowFastaFile, windowCallFile, mode='validation'), 
                    nb_val_samples=nb_validationSamples)