# using the basset architecture to predict enhancer and promoter sites
# basset arch from anna
import numpy as np 
import pdb
from Bio import SeqIO
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from os.path import join, exists
from os import makedirs
from random import shuffle
import sys 

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
    # print('maxEntries ' + str(maxEntries))

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
            numGenerated=0
        batchSequences = np.array([str(fastaGen.next().seq.upper()) for i in range(batchSize)])
        batchOneHot = one_hot_encode(batchSequences)
        batchOneHot = np.transpose(batchOneHot, (0,2,1,3))
        batchCalls = np.array([np.fromstring(cf.readline(), sep='\t',dtype=bool) for i in range(batchSize)])
        
        # take a subset of calls if desired
        if subsetFeatureList is not None:
            batchCalls = batchCalls[:, subsetFeatureList]

        numGenerated += batchSize
        yield tuple([batchOneHot, batchCalls])

# counts number of Ns in window.
# writes output as a newline separated list to outTxt
# otherwise returns count of Ns per window
def count_N_in_window(windowFastaFile, outTxt=None):
    # windowFastaFile = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
    ff = open(windowFastaFile, 'rU')
    fastaGen = SeqIO.parse(ff, "fasta")
    windowN = np.array([])
    while True:
        try:
            seq = fastaGen.next()
            numN = seq.seq.upper().count('N')
            windowN = np.append(windowN, numN)
        except StopIteration:
            break
    if outTxt is not None:
        np.savetxt(outTxt, windowN, delimiter='\n')
    else:
        return(windoN)


# split up the data into train, test, validation sets
# outFiles is a list of 6 files: fasta train, fasta valid, fasta test
# calls train, calls valid, calls test
# if balanced, outputs equal numbers of positive and negatives 
# with this determined by the column balancedColumn, typically the first
def split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles,
 validFraction=0.2, testFraction=0.2, batchSize=10000, balanced=False, balancedColumn=0):
    # get number of fasta sequences in the file
    totalEntries = sum(1 for line in open(inputFasta)) /2
    processedEntries = 0
    totalBatches = int(np.ceil(totalEntries/batchSize))
    print('total batches: ' + str(totalBatches))
    validEntries = int(np.ceil(batchSize * validFraction ))
    testEntries = int(np.ceil(batchSize * testFraction))
    trainEntries = batchSize - validEntries - testEntries

    # open files, read in batches
    with open(inputFasta, 'rU') as ff, open(inputCalls, 'rU') as cf:
        fastaGen = SeqIO.parse(ff, "fasta")
        # read first header line of calls
        cfHeader = cf.readline()

        outFileHandles = [open(i, 'w') for i in outFiles]

        # write header to output Call files
        [i.write(cfHeader) for i in outFileHandles[3:]]
        while processedEntries + batchSize < totalEntries:
            print('Writing batch ' +str(processedEntries/batchSize) + '/' + str(totalBatches))
            # read in a set     
            batchRecords = [fastaGen.next() for i in range(batchSize)]
            # do same thing for calls
            batchCalls = np.array([np.fromstring(cf.readline(), dtype=int, sep='\t') for i in range(batchSize)])
            
            # can do balancing in two ways. From a single task, or
            # counting negatives that are only negative across all classes
            # go with single column for now
            if balanced:
                useCol = batchCalls[:,balancedColumn]
                rowSum = np.sum(batchCalls, axis=1)
                negInd = np.where(useCol==0)[0]
                posInd = np.where(useCol>0)[0]
                npos = len(posInd)
                nneg = len(negInd)
                print('npos: ' + str(npos), ' nneg: ' + str(nneg))
                if nneg > npos:
                    negInd = negInd[:npos]
                elif npos > nneg:
                    posInd = posInd[:nneg]
                finalInd = np.append(np.array(posInd), np.array(negInd))
            else: 
                finalInd = np.array([i for i in range(batchSize)])
            newBatchSize = len(finalInd)
            # print(newBatchSize)

            # shuffle before splitting
            shuffle(finalInd)
            # number of entries in each set
            validEntries = int(np.ceil(newBatchSize * validFraction))
            testEntries = int(np.ceil(newBatchSize * testFraction))
            trainEntries = newBatchSize - validEntries - testEntries
            # get index for each set
            trainIdx = np.sort(finalInd[0:trainEntries])
            validIdx = np.sort(finalInd[trainEntries:trainEntries+validEntries])
            testIdx = np.sort(finalInd[trainEntries+validEntries:])
            # split up records
            trainRecords = [batchRecords[i] for i in trainIdx]
            validRecords = [batchRecords[i] for i in validIdx]
            testRecords = [batchRecords[i] for i in testIdx]
            # write out
            if trainEntries > 1:
                SeqIO.write(trainRecords, outFileHandles[0], 'fasta')
                trainBatch = '\n'.join(['\t'.join([str(a) for a in b]) for b in batchCalls[trainIdx,]]) + '\n'
                outFileHandles[3].write(trainBatch)
            if validEntries > 1:
                SeqIO.write(validRecords, outFileHandles[1], 'fasta')
                validBatch = '\n'.join(['\t'.join([str(a) for a in b]) for b in batchCalls[validIdx,]]) + '\n'
                outFileHandles[4].write(validBatch)
            if testEntries > 1:
                SeqIO.write(testRecords, outFileHandles[2], 'fasta')
                testBatch = '\n'.join(['\t'.join([str(a) for a in b]) for b in batchCalls[testIdx,]]) + '\n'
                outFileHandles[5].write(testBatch)
            processedEntries += batchSize

        [i.close() for i in outFileHandles]

# code to split up inputs
# inputFasta = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
# inputCalls = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls/all.calls'
# outFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
# outFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa','calls_train.txt','calls_valid.txt','calls_test.txt']
# outFiles = [join(outFileBase, i) for i in outFileAdd]
# split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles)
# for making a balanced dataset
# for eash task 
# tasks = [i for i in range(32)]
# for task in tasks:
#     inputFasta = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
#     inputCalls = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls/all.calls'
#     outFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test_balanced/task_' + str(task)
#     if not exists(outFileBase):
#         print('making ' + outFileBase)
#         makedirs(outFileBase)

#     outFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa','calls_train.txt','calls_valid.txt','calls_test.txt']
#     outFiles = [join(outFileBase, i) for i in outFileAdd]
#     split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles, balanced=True, balancedColumn=task)


# metrics from old Keras code
# see https://github.com/keras-team/keras/blob/53e541f7bf55de036f4f5641bd2947b96dd8c4c3/keras/metrics.py
def precision(y_true, y_pred):
    from keras import backend as K
    """Precision metric.
    Only computes a batch-wise average of precision.
    Computes the precision, a metric for multi-label classification of
    how many selected items are relevant.
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)), axis=-1)
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)), axis=-1)
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision


def recall(y_true, y_pred):
    from keras import backend as K
    """Recall metric.
    Only computes a batch-wise average of recall.
    Computes the recall, a metric for multi-label classification of
    how many relevant items are selected.
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)), axis=-1)
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)), axis=-1)
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def positive_accuracy(y_true, y_pred):
    from keras import backend as K
    """
    accuracy only in positive class
    true positives / possible positives 
    """
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    pa = true_positives / (possible_positives + K.epsilon())
    return pa

def negative_accuracy(y_true, y_pred):
    from keras import backend as K
    """
    accuracy only in negative class
    true negatives / possible negatives 
    """
    predicted_negatives = K.sum(K.round(1-y_pred))

    true_negatives = K.sum(K.round((1 - y_true) * (1 - y_pred)))
    possible_negatives = K.sum(K.round((1- y_true)))
    na = true_negatives / (possible_negatives + K.epsilon())
    return na

def predicted_positives(y_true, y_pred):
    from keras import backend as K
    return K.sum(K.round(K.clip(y_pred, 0, 1)))
def predicted_negatives(y_true, y_pred):
    from keras import backend as K
    return K.sum(K.round(1-y_pred))
def possible_positives(y_true, y_pred):
    from keras import backend as K
    return K.sum(K.round(y_true))
def possible_negatives(y_true, y_pred):
    from keras import backend as K
    return K.sum(K.round(1- y_true))

def true_negatives(y_true, y_pred):
    from keras import backend as K
    return K.sum(K.round(K.clip(y_pred, 0, 1)), axis=-1)

def get_weighted_binary_crossentropy(w0_weights, w1_weights):
    import numpy as np
    from keras import backend as K
    # Compute the task-weighted cross-entropy loss, where every task is weighted by 1 - (fraction of non-ambiguous examples that are positive)
    # In addition, weight everything with label -1 to 0
    w0_weights=np.array(w0_weights);
    w1_weights=np.array(w1_weights);
    thresh=-0.5
    def weighted_binary_crossentropy(y_true,y_pred):
        weightsPerTaskRep = y_true*w1_weights[None,:] + (1-y_true)*w0_weights[None,:]
        nonAmbig = K.cast((y_true > -0.5),'float32')
        nonAmbigTimesWeightsPerTask = nonAmbig * weightsPerTaskRep
        return K.mean(K.binary_crossentropy(y_pred, y_true)*nonAmbigTimesWeightsPerTask, axis=-1);
    return weighted_binary_crossentropy; 


# nfeature is the number of output features
def getModelGivenModelOptionsAndWeightInits(w0_file,w1_file,init_weights, taskSubset=None):
    import numpy as np
    np.random.seed(1234)
    import keras
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

    # taskSubset can be used to select a list of tasks for the classifier
    # if none, then do all, 32 in this case
    if taskSubset is None:
        taskSubset = np.array([i for i in range(32)])
        nTask = len(taskSubset)
    elif type(taskSubset) == int:
        nTask = 1
    else:
        nTask = len(taskSubset)

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

    model.add(Dense(nTask))
    model.add(Activation("sigmoid"))
    metrics=["accuracy",precision,recall, positive_accuracy,
            negative_accuracy, predicted_positives, predicted_negatives,
            possible_positives, possible_negatives]
    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    if (w0_file!=None and w1_file!=None):
        w0=[float(i) for i in open(w0_file,'r').read().strip().split('\n')]
        w1=[float(i) for i in open(w1_file,'r').read().strip().split('\n')]
        if nTask==1:
            w0=[w0[taskSubset]]
            w1=[w1[taskSubset]]
        else:
            w0=[w0[i] for i in taskSubset]
            w1=[w1[i] for i in taskSubset]
        loss=get_weighted_binary_crossentropy(w0_weights=w0,w1_weights=w1)
        model.compile(optimizer=adam,loss=loss,
            metrics=metrics)
    else:
        model.compile(optimizer=adam,loss="binary_crossentropy",
            metrics=metrics)
        # pdb.set_trace()
    return model

if __name__ == '__main__':
    main()
