# using the basset architecture to predict enhancer and promoter sites
# basset arch from anna
import numpy as np 
import pdb
from Bio import SeqIO
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

# one hot encode from dragonn repo
def one_hot_encode(sequences):
    sequence_length = len(sequences[0])
    integer_type = np.int8 if sys.version_info[
        0] == 2 else np.int32  # depends on Python version
    integer_array = LabelEncoder().fit(np.array(('ACGTN',)).view(integer_type)).transform(
        sequences.view(integer_type)).reshape(len(sequences), sequence_length)
    one_hot_encoding = OneHotEncoder(
        sparse=False, n_values=5, dtype=integer_type).fit_transform(integer_array)

    return one_hot_encoding.reshape(
        len(sequences), 1, sequence_length, 5).swapaxes(2, 3)[:, :, [0, 1, 2, 4], :]


# data generator for fasta sequences and the corresponding 
# cell type enhancer and promoter calls
# subsetFeatureList can be a np array of indices to use from the call file
def data_generator(windowFastaFile, windowCallFile, batchSize=256, callHeader=True, subsetFeatureList=None):
    numGenerated=0
    # get number of fasta sequences in the file
    totalEntries = sum(1 for line in open(windowCallFile)) - int(callHeader)
    maxEntries = int(totalEntries / batchSize) * batchSize

    ff = open(windowFastaFile, 'rU')
    cf = open(windowCallFile, 'rU')
    fastaGen = SeqIO.parse(ff, "fasta")
    # read first header line of calls
    if callHeader:
        trash=cf.readline()
        # callHeader = np.fromstring(cf.readline(), sep='\t', dtype=str)
    while True:
        if (numGenerated >= maxEntries):
            # reset file handles
            ff.close()
            cf.close()
            ff = open(windowFastaFile, 'rU')
            cf = open(windowCallFile, 'rU')
            fastaGen = SeqIO.parse(ff, "fasta")
            # read first header line of calls
            if callHeader:
                trash=cf.readline()

        batchSequences = np.array([str(fastaGen.next().seq.upper()) for i in range(batchSize)])
        batchOneHot = one_hot_encode(batchSequences)
        batchOneHot = np.transpose(batchOneHot, (0,2,1,3))
        batchCalls = np.array([np.fromstring(cf.readline(), sep='\t',dtype=bool) for i in range(batchSize)])
        
        # take a subset of calls if desired
        if subsetFeatureList is not None:
            batchCalls = batchCalls[:, subsetFeatureList]

        numGenerated += batchSize
        yield tuple([batchOneHot, batchCalls])


# split up the data into train, test, validation sets
# outFiles is a list of 6 files: fasta train, fasta valid, fasta test
# calls train, calls valid, calls test
def split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles,
 validFraction=0.2, testFraction=0.2, batchSize=10000):
    # get number of fasta sequences in the file
    totalEntries = sum(1 for line in open(inputFasta)) /2
    writtenEntries = 0
    totalBatches = int(np.ceil(totalEntries/batchSize))
    validEntries = int(np.ceil(batchSize * validFraction ))
    testEntries = int(np.ceil(batchSize * testFraction))
    trainEntries = batchSize - validEntries - testEntries

    # open files, read in batches
    with open(windowFastaFile, 'rU') as ff, open(windowCallFile, 'rU') as cf:
        fastaGen = SeqIO.parse(ff, "fasta")
        # read first header line of calls
        cfHeader = cf.readline()

        outFileHandles = [open(i, 'w') for i in outFiles]

        # write header to output Call files
        [i.write(cfHeader) for i in outFileHandles[3:]]
        while writtenEntries + batchSize < totalEntries:
            print('Writing batch ' +str(writtenEntries/batchSize) + '/' + str(totalBatches))
            # read in a set     
            batchRecords = [fastaGen.next() for i in range(batchSize)]
            # create random idx
            randIdx = [i for i in range(batchSize)]
            random.shuffle(randIdx)
            trainIdx = np.sort(randIdx[0:trainEntries])
            validIdx = np.sort(randIdx[trainEntries:trainEntries+validEntries])
            testIdx = np.sort(randIdx[trainEntries+validEntries:])
            # split up records
            trainRecords = [batchRecords[i] for i in trainIdx]
            validRecords = [batchRecords[i] for i in validIdx]
            testRecords = [batchRecords[i] for i in testIdx]
            # write out
            SeqIO.write(trainRecords, outFileHandles[0], 'fasta')
            SeqIO.write(validRecords, outFileHandles[1], 'fasta')
            SeqIO.write(testRecords, outFileHandles[2], 'fasta')
            # do same thing for calls
            batchCalls = np.array([cf.readline() for i in range(batchSize)])
            trainBatch = batchCalls[trainIdx]
            validBatch = batchCalls[validIdx]
            testBatch = batchCalls[testIdx]
            # write
            outFileHandles[3].write(trainBatch)
            outFileHandles[4].write(validBatch)
            outFileHandles[5].write(testBatch)
            writtenEntries += batchSize

        # one final case for the last few
        # print('Writing last batch')
        # batchSize = (totalEntries -writtenEntries)
        # # read in a set     
        # batchRecords = [fastaGen.next() for i in range(batchSize)]
        # # create random idx
        # randIdx = [i for i in range(batchSize)]
        # random.shuffle(randIdx)
        # trainIdx = np.sort(randIdx[0:trainEntries])
        # validIdx = np.sort(randIdx[trainEntries:trainEntries+validEntries])
        # testIdx = np.sort(randIdx[trainEntries+validEntries:])
        # # split up records
        # trainRecords = [batchRecords[i] for i in trainIdx]
        # validRecords = [batchRecords[i] for i in validIdx]
        # testRecords = [batchRecords[i] for i in testIdx]
        # # write out
        # SeqIO.write(trainRecords, outFileHandles[0], 'fasta')
        # SeqIO.write(validRecords, outFileHandles[1], 'fasta')
        # SeqIO.write(testRecords, outFileHandles[2], 'fasta')
        # # do same thing for calls
        # batchCalls = np.array([cf.readline() for i in range(batchSize)])
        # trainBatch = batchCalls[trainIdx]
        # validBatch = batchCalls[validIdx]
        # testBatch = batchCalls[testIdx]
        # # write
        # outFileHandles[3].write(trainBatch)
        # outFileHandles[4].write(validBatch)
        # outFileHandles[5].write(testBatch)
        # writtenEntries += batchSize

        [i.close() for i in outFileHandles]

# metrics from old Keras code
# see https://github.com/keras-team/keras/blob/53e541f7bf55de036f4f5641bd2947b96dd8c4c3/keras/metrics.py
def precision(y_true, y_pred):
    """Precision metric.
    Only computes a batch-wise average of precision.
    Computes the precision, a metric for multi-label classification of
    how many selected items are relevant.
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def recall(y_true, y_pred):
    """Recall metric.
    Only computes a batch-wise average of recall.
    Computes the recall, a metric for multi-label classification of
    how many relevant items are selected.
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

# code to split up inputs
# inputFasta = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
# inputCalls = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls/all.calls'
# outFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
# outFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa','calls_train.txt','calls_valid.txt','calls_test.txt']
# outFiles = [join(outFileBase, i) for i in outFileAdd]
# split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles)
# nfeature is the number of output features
def getModelGivenModelOptionsAndWeightInits(w0_file,w1_file,init_weights, nFeature=32):
    import numpy as np
    np.random.seed(1234)
    import keras;
    from keras.models import Sequential
    from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
    from keras.layers.convolutional import Conv2D, MaxPooling2D
    from keras.optimizers import Adadelta, SGD, RMSprop;
    # from keras.losses import *;
    from keras.constraints import maxnorm;
    from keras.layers.normalization import BatchNormalization
    from keras.regularizers import l1, l2    
    from keras import backend as K
    from keras import metrics

    model=Sequential()
    model.add(Conv2D(filters=300,kernel_size=(1,19),input_shape=(4,1,2000)))

    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1,4)))

    model.add(Conv2D(filters=200,kernel_size=(1,11)))
    model.add(BatchNormalization(axis=1))
    model.add(Activation('relu'))
    model.add(MaxPooling2D(pool_size=(1,4)))

    model.add(Conv2D(filters=200,kernel_size=(1,7)))
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

    model.add(Dense(nFeature))
    model.add(Activation("sigmoid"))

    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    if (w0_file!=None and w1_file!=None):
        w0=[float(i) for i in open(w0_file,'r').read().strip().split('\n')]
        w1=[float(i) for i in open(w1_file,'r').read().strip().split('\n')]
        loss=get_weighted_binary_crossentropy(w0_weights=w0,w1_weights=w1)
        model.compile(optimizer=adam,loss=loss,
            metrics=["accuracy"])
    else:
        model.compile(optimizer=adam,loss="binary_crossentropy",
            metrics=["accuracy", precision,recall])
            # metrics=["accuracy"])
        # pdb.set_trace()
    return model

'''
# and do same prediction as before
from os.path import join, exists
from os import makedirs, chdir
import numpy as np
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
'''


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

'''
model=getModelGivenModelOptionsAndWeightInits(None, None, None)

inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa','calls_train.txt','calls_valid.txt','calls_test.txt']
inFiles = [join(inFileBase, i) for i in inFileAdd]

nb_trainSamples = int((sum(1 for line in open(inFiles[0])) /2) * 0.98)
nb_validationSamples = int((sum(1 for line in open(inFiles[1])) /2) * 0.98)

# or for testing 
nb_trainSamples = 100000
nb_validationSamples = 100000

model.fit_generator(generator=data_generator(inFiles[0], inFiles[3], batchSize=512),
                    steps_per_epoch =  nb_trainSamples, epochs=100, show_accuracy=True,
                    validation_data = data_generator(inFiles[1], inFiles[4], batchSize=512), 
                    validation_steps=validSteps)

'''


# only predicting a few features, should be easier
# useFeatures = np.array([0,2,4,6])
useFeatures = np.array([i for i in range(32)])
model=getModelGivenModelOptionsAndWeightInits(None, None, None, nFeature=len(useFeatures))
from os.path import join
import numpy as np
import sys
inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa',
             'calls_train.txt','calls_valid.txt','calls_test.txt']
# inFileAdd = ['short_fasta_train.fa','short_fasta_valid.fa','short_fasta_test.fa',
#              'short_calls_train.txt','short_calls_valid.txt','short_calls_test.txt']
inFiles = [join(inFileBase, i) for i in inFileAdd]

batchSize=256
nb_trainSamples = int(sum(1 for line in open(inFiles[3])) -1)
nb_validationSamples = int(sum(1 for line in open(inFiles[4])) -1)
trainSteps = np.floor(nb_trainSamples / batchSize)
validSteps = np.floor(nb_trainSamples / batchSize)

# or for testing 
# nb_trainSamples = 10240
# nb_validationSamples = 1024

model.fit_generator(generator=data_generator(inFiles[0], inFiles[3], batchSize=batchSize, subsetFeatureList=useFeatures),
                    steps_per_epoch =  trainSteps, epochs=10, 
                    validation_data = data_generator(inFiles[1], inFiles[4], batchSize=batchSize, subsetFeatureList=useFeatures), 
                    validation_steps=validSteps)


