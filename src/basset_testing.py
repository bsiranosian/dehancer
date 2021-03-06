# using the basset architecture to predict enhancer and promoter sites
# basset arch from anna
import numpy as np 
import pdb
from Bio import SeqIO
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from os.path import join
from random import shuffle

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
            
            # now do balancing by counting negatives that are only negative across
            # all classes
            if balanced:
                rowSum = np.sum(batchCalls, axis=1)
                negInd = np.where(rowSum==0)[0]
                posInd = np.where(rowSum>0)[0]
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

            validEntries = int(np.ceil(newBatchSize * validFraction ))
            testEntries = int(np.ceil(newBatchSize * testFraction))
            trainEntries = newBatchSize - validEntries - testEntries


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
# inputFasta = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/fasta/genome_tile_windows.fa'
# inputCalls = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/window_calls/all.calls'
# outFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test_balanced'
# outFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa','calls_train.txt','calls_valid.txt','calls_test.txt']
# outFiles = [join(outFileBase, i) for i in outFileAdd]
# split_fasta_calls_train_validate_test(inputFasta, inputCalls, outFiles, balanced=True)


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
def getModelGivenModelOptionsAndWeightInits(w0_file,w1_file,init_weights, nFeature=32):
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
    metrics=["accuracy",precision,recall, positive_accuracy,
            negative_accuracy, predicted_positives, predicted_negatives,
            possible_positives, possible_negatives]
    adam = keras.optimizers.Adam(lr=0.001, beta_1=0.9, beta_2=0.999, epsilon=1e-08)
    if (w0_file!=None and w1_file!=None):
        w0=[float(i) for i in open(w0_file,'r').read().strip().split('\n')]
        w1=[float(i) for i in open(w1_file,'r').read().strip().split('\n')]
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

def main():
    from os.path import join, exists
    from os import makedirs
    import numpy as np
    import sys

    balanced =True
    if not balanced:
        allFeatures=True
        if allFeatures:
            # enhancers and promoters across 16 cell lines
            useFeatures = np.array([i for i in range(32)])
            # positive and negative call weights, for each cell line
            positiveCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/positive_call_weights.txt'
            negativeCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/negative_call_weights.txt'

            trainEpochs = 10
            saveDir = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/allcell_enh_pro/'
            saveName = 'basset_model_trained_' + str(trainEpochs) + '.hdf5'

        else: 
            # only predicting a few features, should be easier
            # Enhancers in Dnd41 only
            useFeatures = np.array([0])
            # positive and negative call weights, for each cell line
            positiveCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/Dnd41_positive_call_weights.txt'
            negativeCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/Dnd41_negative_call_weights.txt'

            trainEpochs = 10
            saveDir = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/Dnd41_enh/'
            saveName = 'basset_model_trained_' + str(trainEpochs) + '.hdf5'


        # create save dir if it doesnt exist
        if not exists(saveDir):
            print('making ' + saveDir)
            makedirs(saveDir)
        # create model
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
        inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa',
                     'calls_train.txt','calls_valid.txt','calls_test.txt']
        # inFileAdd = ['short_fasta_train.fa','short_fasta_valid.fa','short_fasta_test.fa',
        #              'short_calls_train.txt','short_calls_valid.txt','short_calls_test.txt']
        inFiles = [join(inFileBase, i) for i in inFileAdd]


        model=getModelGivenModelOptionsAndWeightInits(
            positiveCallWeightsFile, negativeCallWeightsFile,
            None, nFeature=len(useFeatures))

        batchSize=256
        nb_trainSamples = int(sum(1 for line in open(inFiles[3])) -1)
        nb_validationSamples = int(sum(1 for line in open(inFiles[4])) -1)
        trainSteps = np.floor(nb_trainSamples / batchSize)
        validSteps = np.floor(nb_trainSamples / batchSize)

        trainSteps=10
        validSteps=10
        trainGenerator = data_generator(inFiles[0], inFiles[3],
            batchSize=batchSize, subsetFeatureList=useFeatures)
        validGenerator = data_generator(inFiles[1], inFiles[4], 
            batchSize=batchSize, subsetFeatureList=useFeatures)

        # train model
        history  = model.fit_generator(generator=trainGenerator,
            steps_per_epoch=trainSteps, epochs=trainEpochs, 
            validation_data=validGenerator, validation_steps=validSteps)

        model.save(join(saveDir, saveName))

    elif balanced:
        # balanced dataset
        allFeatures=True
        if allFeatures:
            # enhancers and promoters across 16 cell lines
            useFeatures = np.array([i for i in range(32)])
            trainEpochs = 100
            saveDir = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/allcell_enh_pro_balanced/'
            saveDirFigures = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/allcell_enh_pro_balanced/figures'
            saveName = 'basset_model_trained_' + str(trainEpochs) + '.hdf5'
            modelFile = join(saveDir, saveName)
        else: 
            # only predicting a few features, should be easier
            # Enhancers in Dnd41 only
            useFeatures = np.array([0])
            trainEpochs = 100
            saveDir = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/Dnd41_enh_balanced/'
            saveDirFigures = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/Dnd41_enh_balanced/figures'
            saveName = 'basset_model_trained_' + str(trainEpochs) + '.hdf5'
            modelFile = join(saveDir, saveName)


        # create save dir if it doesnt exist
        if not exists(saveDir):
            print('making ' + saveDir)
            makedirs(saveDir)
        if not exists(saveDirFigures):
            print('making ' + saveDirFigures)
            makedirs(saveDirFigures)
        # create model
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test_balanced'
        inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa',
                     'calls_train.txt','calls_valid.txt','calls_test.txt']
        # inFileAdd = ['short_fasta_train.fa','short_fasta_valid.fa','short_fasta_test.fa',
        #              'short_calls_train.txt','short_calls_valid.txt','short_calls_test.txt']
        inFiles = [join(inFileBase, i) for i in inFileAdd]



        batchSize=256
        nb_trainSamples = int(sum(1 for line in open(inFiles[3])) -1)
        nb_validationSamples = int(sum(1 for line in open(inFiles[4])) -1)
        nb_testSamples = int(sum(1 for line in open(inFiles[5])) -1)
        trainSteps = np.floor(nb_trainSamples / batchSize)
        validSteps = np.floor(nb_validationSamples / batchSize)
        testSteps = np.floor(nb_testSamples / batchSize)

        trainGenerator = data_generator(inFiles[0], inFiles[3],
            batchSize=batchSize, subsetFeatureList=useFeatures)
        validGenerator = data_generator(inFiles[1], inFiles[4], 
            batchSize=batchSize, subsetFeatureList=useFeatures)
        testGenerator = data_generator(inFiles[2], inFiles[5], 
            batchSize=batchSize, subsetFeatureList=useFeatures)

        trainModel=True
        # train model if desired
        if trainModel:
            model=getModelGivenModelOptionsAndWeightInits(
                None, None,
                None, nFeature=len(useFeatures))
            history  = model.fit_generator(generator=trainGenerator,
                steps_per_epoch=trainSteps, epochs=trainEpochs, 
                validation_data=validGenerator, validation_steps=validSteps)
            # save model
            model.save(modelFile)

        else:
            # load model if already trained
            custom_objects = {'precision': precision,
            'recall': recall,
            'positive_accuracy': positive_accuracy,
            'negative_accuracy': negative_accuracy,
            'predicted_positives': predicted_positives,
            'predicted_negatives': predicted_negatives,
            'possible_positives': possible_positives,
            'possible_negatives': possible_negatives}
            model = keras.models.load_model(modelFile, custom_objects=custom_objects)

        # evaluation on a test dataset    
        testPredict  = model.predict_generator(generator=testGenerator, steps=testSteps)
        testTrue = np.loadtxt(inFiles[5], dtype=bool, skiprows=1, delimiter='\t')[0:int(testSteps * batchSize),]
        with open(inFiles[5], 'rU') as inf:
            testHeader = inf.readline().strip().split('\t')

        # get precision recall curve
        import matplotlib.pyplot as plt
        from sklearn.metrics import precision_recall_curve, average_precision_score

        ## which task to test on
        for task in range(32):  
            taskName = testHeader[task]
            prcSaveName = join(saveDirFigures, taskName + '_prc.png')

            pr, re, _ =  precision_recall_curve(testTrue[:,task], testPredict[:,task])
            average_precision = average_precision_score(testTrue[:,task], testPredict[:,task], average="micro")

            plt.step(re, pr, color='b', alpha=0.2,
                     where='post')
            plt.fill_between(re, pr, step='post', alpha=0.2,
                             color='b')

            plt.xlabel('Recall')
            plt.ylabel('Precision')
            plt.ylim([0.0, 1.05])
            plt.xlim([0.0, 1.0])
            plt.title(taskName + ' PRC | AP={0:0.2f}'.format(
                      average_precision))
            # plt.show(block=False)
            plt.savefig(prcSaveName)
            plt.close()