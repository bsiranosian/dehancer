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
from os import makedirs
import numpy as np
# single task across all tasks, unbalanced dataset
# input files
positiveCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/positive_call_weights.txt'
negativeCallWeightsFile='/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/negative_call_weights.txt'

#output dirs
baseDir = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/trained_models/single_task_models/'

# info about different tasks
tasks = np.array([i for i in range(32)])
with open('/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test/calls_train.txt', 'rU') as inf:
    taskHeader = inf.readline().strip().split('\t')
# model settings
trainEpochs = 25
batchSize = 256
weighted=False
balanced=True
for task in tasks:
    task=int(task)
    taskName = taskHeader[task]
    print('Starting task ' + taskName)

    if balanced:
        saveDir = join(baseDir, taskName + '_balanced')
    else:
        saveDir = join(baseDir, taskName + '_unbalanced')

    modelSaveName = 'basset_model_trained_' + str(trainEpochs) + '.hdf5'

    # create save dir if it doesnt exist
    if not exists(saveDir):
        print('making ' + saveDir)
        makedirs(saveDir)

    # input files
    if balanced:
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test_balanced/task_' + str(task)
    else:
        inFileBase = '/users/bsiranos/analysis/enhancer_conservation/encode_data/broad/ml_train_valid_test'
    inFileAdd = ['fasta_train.fa','fasta_valid.fa','fasta_test.fa',
                 'calls_train.txt','calls_valid.txt','calls_test.txt']
    inFiles = [join(inFileBase, i) for i in inFileAdd]
    # init model
    if weighted and not balanced:
        model=getModelGivenModelOptionsAndWeightInits(
            positiveCallWeightsFile, negativeCallWeightsFile,
            None, taskSubset=task)
    else:
        model=getModelGivenModelOptionsAndWeightInits(
            None, None, None, taskSubset=task)

    nb_trainSamples = int(sum(1 for line in open(inFiles[3])) -1)
    nb_validationSamples = int(sum(1 for line in open(inFiles[4])) -1)
    nb_testSamples = int(sum(1 for line in open(inFiles[5])) -1)
    trainSteps = np.floor(nb_trainSamples / batchSize)
    validSteps = np.floor(nb_validationSamples / batchSize)
    testSteps = np.floor(nb_testSamples / batchSize)

    # trainSteps=10
    # validSteps=10
    trainGenerator = data_generator(inFiles[0], inFiles[3],
        batchSize=batchSize, subsetFeatureList=task)
    validGenerator = data_generator(inFiles[1], inFiles[4], 
        batchSize=batchSize, subsetFeatureList=task)
    testGenerator = data_generator(inFiles[2], inFiles[5], 
        batchSize=batchSize, subsetFeatureList=task)

    # train model
    history  = model.fit_generator(generator=trainGenerator,
        steps_per_epoch=trainSteps, epochs=trainEpochs, 
        validation_data=validGenerator, validation_steps=validSteps)

    model.save(join(saveDir, modelSaveName))

    # evaluation on a test dataset    
    testPredict  = model.predict_generator(generator=testGenerator, steps=testSteps)
    testTrue = np.loadtxt(inFiles[5], dtype=bool, skiprows=1, delimiter='\t')[0:int(testSteps * batchSize), task]

    # get precision recall curve
    import matplotlib.pyplot as plt
    from sklearn.metrics import precision_recall_curve, average_precision_score

    # prediction evaluation

    ## which task to test on
    # for task in range(32):  


    pr, re, _ =  precision_recall_curve(testTrue, testPredict)
    average_precision = average_precision_score(testTrue, testPredict, average="micro")

    prcSaveName = join(saveDir, taskName + '_prc.png')
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


'''
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

# prediction evaluation

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
'''
