### https://github.com/laserwave/plsa
import numpy as np
from numpy import zeros, int8, log
# from pylab import random
import pickle
import random
import sys
import jieba
import re
import time
import codecs
import matplotlib.pyplot as plt
# import ipdb


def initializeParameters():
    for i in range(0, N):
        normalization = sum(lamda[i, :])
        for j in range(0, K):
            lamda[i, j] /= normalization;

    for i in range(0, K):
        normalization = sum(theta[i, :])
        for j in range(0, M):
            theta[i, j] /= normalization;

    
def preprocessing(datasetFilePath, stopwordsFilePath):
    
    # read the stopwords file
    file = codecs.open(stopwordsFilePath, 'r', 'utf-8')
    stopwords = [line.strip() for line in file] 
    file.close()
    
    # read the documents
    file = codecs.open(datasetFilePath, 'r', 'utf-8')
    # file = codecs.open(datasetFilePath, 'r', encoding="latin-1")
    documents = [document.strip() for document in file] 
    file.close()

    # number of documents
    N = len(documents)

    wordCounts = [];
    word2id = {}
    id2word = {}
    currentId = 0;
    # generate the word2id and id2word maps and count the number of times of words showing up in documents
    for document in documents:
        segList = jieba.cut(document)
        wordCount = {}
        for word in segList:
            word = word.lower().strip()
            if len(word) > 1 and not re.search('[0-9]', word) and word not in stopwords:               
                if word not in word2id.keys():
                    word2id[word] = currentId;
                    id2word[currentId] = word;
                    currentId += 1;
                if word in wordCount:
                    wordCount[word] += 1
                else:
                    wordCount[word] = 1
        wordCounts.append(wordCount);
    
    # length of dictionary
    M = len(word2id)  

    # generate the document-word matrix
    X = zeros([N, M], int8)
    for word in word2id.keys():
        j = word2id[word]
        for i in range(0, N):
            if word in wordCounts[i]:
                X[i, j] = wordCounts[i][word];    

    return N, M, word2id, id2word, X

def EStep():
    for i in range(0, N):
        for j in range(0, M):
            denominator = 0;
            for k in range(0, K):
                p[i, j, k] = theta[k, j] * lamda[i, k];
                denominator += p[i, j, k];
            if denominator == 0:
                for k in range(0, K):
                    p[i, j, k] = 0;
            else:
                for k in range(0, K):
                    p[i, j, k] /= denominator;



def EStep_incremental(index):
    for i in range(0, N):
        for j in range(0, M):
            denominator = 0;
            for k in range(0, K):
                if i in index:
                    p[i, j, k] = theta[k, j] * lamda[i, k];
                else: 
                    p[i, j, k] = oldp[i, j, k]
                denominator += p[i, j, k];
            if denominator == 0:
                for k in range(0, K):
                    p[i, j, k] = 0;
            else:
                for k in range(0, K):
                    p[i, j, k] /= denominator;


def MStep():
    # update theta
    for k in range(0, K):
        denominator = 0
        for j in range(0, M):
            theta[k, j] = 0
            for i in range(0, N):
                theta[k, j] += X[i, j] * p[i, j, k]
            denominator += theta[k, j]
        if denominator == 0:
            for j in range(0, M): 
                theta[k, j] = 1.0 / M
        else:
            for j in range(0, M):
                theta[k, j] /= denominator
        
    # update lamda
    for i in range(0, N):
        for k in range(0, K):
            lamda[i, k] = 0
            denominator = 0
            for j in range(0, M):
                lamda[i, k] += X[i, j] * p[i, j, k]
                # lamda[i, k] += 1
                denominator += X[i, j];
            if denominator == 0:
                lamda[i, k] = 1.0 / K
            else:
                lamda[i, k] /= denominator

def MStep_online(index):
    # update theta
    for k in range(0, K):
        denominator = 0
        for j in range(0, M):
            theta[k, j] = 0
            oldsomme = 0
            somme_minibatch = 0
            for i in range(0,N):
                oldsomme += X[i, j] * oldp[i, j, k]
            for i in index:
                somme_minibatch += X[i, j] * p[i, j, k]        
            theta[k, j] += oldsomme + rho[epoch]*(somme_minibatch - oldsomme)
            denominator += theta[k, j]
        if denominator == 0:
            for j in range(0, M):
                theta[k, j] = 1.0 / M
        else:
            for j in range(0, M):
                theta[k, j] /= denominator
        
    # update lamda
    for i in range(0, N):
        for k in range(0, K):
            lamda[i, k] = 0
            denominator = 0
            for j in range(0, M):
                lamda[i, k] += X[i, j] * p[i, j, k]
                denominator += X[i, j];
            if denominator == 0:
                lamda[i, k] = 1.0 / K
            else:
                lamda[i, k] /= denominator

def MStep_onlinevr(index):
    # update theta
    for k in range(0, K):
        denominator = 0
        for j in range(0, M):
            theta[k, j] = 0
            oldsomme = 0
            oldsomme0 = 0
            somme_minibatch = 0
            somme_minibatch0 = 0
            for i in range(0,N):
                oldsomme += X[i, j] * oldp[i, j, k]
                oldsomme0 += X[i, j] * p0[i, j, k]
            for i in index:
                somme_minibatch += X[i, j] * p[i, j, k]
                somme_minibatch0 += X[i, j] * p0[i, j, k]
            theta[k, j] += oldsomme + rho*(somme_minibatch - somme_minibatch0 + oldsomme0 - oldsomme)
            denominator += theta[k, j]
        if denominator == 0:
            for j in range(0, M):
                theta[k, j] = 1.0 / M
        else:
            for j in range(0, M):
                theta[k, j] /= denominator
        
    # update lamda
    for i in range(0, N):
        for k in range(0, K):
            lamda[i, k] = 0
            denominator = 0
            for j in range(0, M):
                lamda[i, k] += X[i, j] * p[i, j, k]
                denominator += X[i, j];
            if denominator == 0:
                lamda[i, k] = 1.0 / K
            else:
                lamda[i, k] /= denominator


def SAGAStep(index_i, index_j):
    # compute new and old individual h_i
    oldh = h
    for i in index_i:
        for j in range(0, M):
            denominator = 0;
            for k in range(0, K):
                h[i, j, k] = theta[k, j] * lamda[i, k];
                denominator += h[i, j, k];
            if denominator == 0:
                for k in range(0, K):
                    h[i, j, k] = 0;
            else:
                for k in range(0, K):
                    h[i, j, k] /= denominator;
    
    # update theta
    for k in range(0, K):
        denominator = 0
        for j in range(0, M):    
            tmp = 0
            for i in index_i:
                tmp += X[i, j] *(h[i, j, k] - oldh[i, j, k])
            Vsomme[k,j] = Hsomme[k,j] + N*tmp
            Ssomme[k,j] = (1-rho_saga)*Ssomme[k,j] + rho_saga*Vsomme[k,j]
            theta[k, j] = Ssomme[k, j]
            denominator += theta[k, j]
        if denominator == 0:
            for j in range(0, M): 
                theta[k, j] = 1.0 / M
        else:
            for j in range(0, M):
                theta[k, j] /= denominator

    # update lamda
    for i in range(0, N):
        for k in range(0, K):
            lamda[i, k] = 0
            denominator = 0
            for j in range(0, M):
                lamda[i, k] += X[i, j] * p[i, j, k]
                denominator += X[i, j];
            if denominator == 0:
                lamda[i, k] = 1.0 / K
            else:
                lamda[i, k] /= denominator

    for k in range(0, K):
        for j in range(0, M):
            for i in index_j:
                Hsomme[k,j] += X[i, j] *(theta[k, j] * lamda[i, k] - listofthetas[i][k, j] * listoflamdas[i, k]);
    #save all parameters per individual
    for i in index_j:
        listofthetas[i] = theta
        for k in range(0, K):
            listoflamdas[i,k] = lamda[i,k]


# def MStepsaga(index_i, index_j):
#     # update theta
#     for i in index_i:
#         for j in range(0, M):
#             denominator = 0;
#             for k in range(0, K):
#                 h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
#                 denominator += h[i, j, k];
#             if denominator == 0:
#                 for k in range(0, K):
#                     h[i, j, k] = 0;
#             else:
#                 for k in range(0, K):
#                     h[i, j, k] /= denominator;
#     for k in range(0, K):
#         denominator = 0
#         for j in range(0, M):
#             theta[k, j] = 0
#             tmp = 0
#             for i in index_i:
#                 tmp +=  X[i, j] * p[i, j, k] - X[i, j] * h[i, j, k]
#             Vsomme = Hsomme[k,j] + N*tmp
#             Ssomme[k,j] = (1- rho_saga)*Ssomme[k,j] + rho_saga*Vsomme
#             theta[k, j] = Ssomme[k,j]
#             denominator += theta[k, j]
#         if denominator == 0:
#             for j in range(0, M): 
#                 theta[k, j] = 1.0 / M
#         else:
#             for j in range(0, M):
#                 theta[k, j] /= denominator
#     # update lamda
#     for i in range(0, N):
#         for k in range(0, K):
#             lamda[i, k] = 0
#             denominator = 0
#             for j in range(0, M):
#                 lamda[i, k] += X[i, j] * p[i, j, k]
#                 # lamda[i, k] += 1
#                 denominator += X[i, j];
#             if denominator == 0:
#                 lamda[i, k] = 1.0 / K
#             else:
#                 lamda[i, k] /= denominator
#     for j in range(0, M):
#         for k in range(0, K):
#             tmph = 0
#             for i in index_j:
#                 h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
#                 tmph +=  (X[i, j] * p[i, j, k] - X[i, j] * h[i, j, k])
#             Hsomme[k,j] += tmph
#     for i in index_j:
#         listofthetas[i] = theta
#         listoflamdas[i,] = lamda[i,]


# calculate the log likelihood
def LogLikelihood():
    loglikelihood = 0
    for i in range(0, N):
        for j in range(0, M):
            tmp = 0
            for k in range(0, K):
                tmp += theta[k, j] * lamda[i, k]
            if tmp > 0:
                loglikelihood += X[i, j] * log(tmp)
    return loglikelihood

# output the params of model and top words of topics to files
def output():
    # document-topic distribution
    file = codecs.open(docTopicDist,'w','utf-8')
    for i in range(0, N):
        tmp = ''
        for j in range(0, K):
            tmp += str(lamda[i, j]) + ' '
        file.write(tmp + '\n')
    file.close()
    
    # topic-word distribution
    file = codecs.open(topicWordDist,'w','utf-8')
    for i in range(0, K):
        tmp = ''
        for j in range(0, M):
            tmp += str(theta[i, j]) + ' '
        file.write(tmp + '\n')
    file.close()
    
    # dictionary
    file = codecs.open(dictionary,'w','utf-8')
    for i in range(0, M):
        file.write(id2word[i] + '\n')
    file.close()
    
    # top words of each topic
    file = codecs.open(topicWords,'w','utf-8')
    for i in range(0, K):
        topicword = []
        ids = theta[i, :].argsort()
        for j in ids:
            topicword.insert(0, id2word[j])
        tmp = ''
        for word in topicword[0:min(topicWordsNum, len(topicword))]:
            tmp += word + ' '
        file.write(tmp + '\n')
    file.close()
    


datasetFilePath = 'dataset1.txt'
# datasetFilePath = 'dataset2.txt'
# datasetFilePath = 'dataset5.txt'

# datasetFilePath = 'dataset100k.txt'
# datasetFilePath = 'dataset10k.txt'
datasetFilePath = 'dataset1k.txt'

stopwordsFilePath = 'stopwords.dic'
K = 10    # number of topic
nb_epochs = 10
maxIteration = 30
threshold = 10.0
topicWordsNum = 10
docTopicDist = 'docTopicDistribution.txt'
topicWordDist = 'topicWordDistribution.txt'
dictionary = 'dictionary.dic'
topicWords = 'topics.txt'
if(len(sys.argv) == 11):
    datasetFilePath = sys.argv[1]
    stopwordsFilePath = sys.argv[2]
    K = int(sys.argv[3])
    maxIteration = int(sys.argv[4])
    threshold = float(sys.argv[5])
    topicWordsNum = int(sys.argv[6])
    docTopicDist = sys.argv[7]
    topicWordDist = sys.argv[8]
    dictionary = sys.argv[9]
    topicWords = sys.argv[10]


# preprocessing
N, M, word2id, id2word, X = preprocessing(datasetFilePath, stopwordsFilePath)
print(N)
print(M)
print(X.shape)
print(len(word2id))

#LIST OF INDICES FOR INCREMENTAL METHODS
seed0 = 333888
indices = [x for x in range(N)]
list_indices_i = []
for epoch in range(0, nb_epochs):
    indices = [x for x in range(N)]
    random.seed(seed0*(epoch+3))
    random.shuffle(indices)
    list_indices_i.append(indices)

list_indices_j = []
for epoch in range(0, nb_epochs):
    indices_j = list_indices_i[epoch]
    random.seed(1999*(epoch+1))
    random.shuffle(indices_j)
    list_indices_j.append(indices_j)

seed0 = 333888
indices = [x for x in range(N)]
list_indices_i = []
for epoch in range(0, nb_epochs):
    indices = [x for x in range(N)]
    random.seed(seed0*(epoch+30))
    random.shuffle(indices)
    list_indices_i.append(indices)


mini_batch_size = 4

#REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])

p = zeros([N, M, K])
initializeParameters()
oldLoglikelihood = 1
newLoglikelihood = 1

### Batch EM
objectiveEM = []
for epoch in range(0, nb_epochs):
    EStep()
    MStep()
    newLoglikelihood = LogLikelihood()
    print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
    # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
    #     break
    objectiveEM.append(newLoglikelihood)
    oldLoglikelihood = newLoglikelihood

# with open('losses/localiemloss', 'wb') as fp: 
#     pickle.dump(objectiveIEM, fp)


#REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])
p = zeros([N, M, K])
initializeParameters()
oldLoglikelihood = 1
newLoglikelihood = 1

### Incremental EM
objectiveIEM = []
for epoch in range(0, nb_epochs):
# for epoch in range(0, 2):
    if epoch == 0:
        EStep()
        MStep()
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveIEM.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood
    else:
        mini_batches = [list_indices_i[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
        for mini_batch in mini_batches:
            oldp = p
            EStep_incremental(mini_batch)
            MStep()
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveIEM.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood



#REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])
p = zeros([N, M, K])
initializeParameters()
oldLoglikelihood = 1
newLoglikelihood = 1

### Online EM
objectiveoEM = []
rho = list(map(lambda x: 3/(x+10), list(range(nb_epochs)))) #STEPSIZE

for epoch in range(0, nb_epochs):
    if epoch == 0:
        EStep()
        MStep()
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveoEM.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood
    else:
        indices = [x for x in range(N)]
        random.shuffle(indices)
        mini_batches = [indices[k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
        for mini_batch in mini_batches:
            oldp = p
            EStep_incremental(mini_batch)
            MStep_online(mini_batch)
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveoEM.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood



#REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])
p = zeros([N, M, K])
initializeParameters()
oldLoglikelihood = 1
newLoglikelihood = 1

### Online EM with VR
objectiveoEM_vr = []
#stepsizes for online EM
rho = 0.1 #STEPSIZE FOR VR

for epoch in range(0, nb_epochs):
    p0 = p
    if epoch == 0:
        EStep()
        MStep()
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveoEM_vr.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood
    else:
        indices = [x for x in range(N)]
        random.shuffle(indices)
        mini_batches = [indices[k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
        for mini_batch in mini_batches:
            oldp = p
            EStep_incremental(mini_batch)
            MStep_onlinevr(mini_batch)
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - old Loglikelihood < threshold):
        #     break
        objectiveoEM_vr.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood



# with open ('lossesjap/emloss', 'rb') as fp:
#     objectiveEM = pickle.load(fp)

# with open('losses/sagaloss', 'wb') as fp: 
#     pickle.dump(objectiveSAGA, fp)

# with open ('lossesjap/sagaloss', 'rb') as fp:
#     objectiveSAGA = pickle.load(fp)
# with open ('lossesjap/iemloss', 'rb') as fp:
#     objectiveIEM = pickle.load(fp)
# with open ('lossesjap/oemloss', 'rb') as fp:
#     objectiveoEM = pickle.load(fp)
# with open ('lossesjap/oemvrloss', 'rb') as fp:
#     objectiveoEM_vr = pickle.load(fp)

# #### PLOTTING #######
# plt.plot(np.arange(nb_epochs), objectiveIEM, label='IEM')
# plt.plot(np.arange(nb_epochs), objectiveEM, label='EM')
# plt.plot(np.arange(nb_epochs), objectiveoEM, label='oEM')
# plt.plot(np.arange(nb_epochs), objectiveoEM_vr, label='oEMVR')
# plt.plot(np.arange(nb_epochs), objectiveSAGA, label='FI-EM')
# leg = plt.legend(fontsize=20,fancybox=True, loc='right')
# leg.get_frame().set_alpha(0.5)
# plt.xlabel('Epoch', fontsize=15)
# plt.ylabel('Objective', fontsize=15)
# plt.show()



with open ('losses1k/emloss', 'rb') as fp:
    objectiveEM_1k = pickle.load(fp)
nb_epochs = len(objectiveEM_1k)
with open ('losses1k/sagaloss', 'rb') as fp:
    objectiveSAGA_1k = pickle.load(fp)
with open ('losses1k/iemloss', 'rb') as fp:
    objectiveIEM_1k = pickle.load(fp)
with open ('losses1k/oemloss', 'rb') as fp:
    objectiveoEM_1k = pickle.load(fp)
with open ('losses1k/oemvrloss', 'rb') as fp:
    objectiveoEM_vr_1k = pickle.load(fp)



#### PLOTTING #######
plt.plot(np.arange(nb_epochs), objectiveIEM, label='IEM')
plt.plot(np.arange(nb_epochs), objectiveEM, label='EM')
plt.plot(np.arange(nb_epochs), objectiveoEM, label='oEM')
plt.plot(np.arange(nb_epochs), objectiveoEM_vr, label='oEMVR')
plt.plot(np.arange(nb_epochs), objectiveSAGA, label='FI-EM')
leg = plt.legend(fontsize=20,fancybox=True, loc='right')
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epoch', fontsize=15)
plt.ylabel('Objective', fontsize=15)
plt.show()



### 10k
with open ('losses10k/emloss', 'rb') as fp:
    objectiveEM = pickle.load(fp)
nb_epochs = len(objectiveEM)
with open ('losses10k/sagaloss', 'rb') as fp:
    objectiveSAGA = pickle.load(fp)
with open ('losses10k/iemloss', 'rb') as fp:
    objectiveIEM = pickle.load(fp)
with open ('losses10k/oemloss', 'rb') as fp:
    objectiveoEM = pickle.load(fp)
with open ('losses10k/oemvrloss', 'rb') as fp:
    objectiveoEM_vr = pickle.load(fp)


#### PLOTTING #######
plt.plot(np.arange(nb_epochs), objectiveIEM, label='IEM')
plt.plot(np.arange(nb_epochs), objectiveEM, label='EM')
plt.plot(np.arange(nb_epochs), objectiveoEM, label='oEM')
plt.plot(np.arange(nb_epochs), objectiveoEM_vr, label='oEMVR')
plt.plot(np.arange(nb_epochs), objectiveSAGA, label='FI-EM')
leg = plt.legend(fontsize=20,fancybox=True, loc='right')
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epoch', fontsize=15)
plt.ylabel('Objective', fontsize=15)
plt.show()





xaxis = np.arange(nb_epochs)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16, 6))
ax = plt.subplot(1, 2, 1)
plt.plot(xaxis, objectiveIEM_1k, label='IEM', marker='^')
plt.plot(xaxis, objectiveEM_1k, label='EM', marker='^')
plt.plot(xaxis, objectiveSAGA_1k, label='FI-EM', marker='^')
plt.plot(xaxis, objectiveoEM_vr_1k, label='SVR-EM', marker='^')
plt.plot(xaxis, objectiveoEM_1k, label='OEM', marker='^')
leg = plt.legend(fontsize=20,fancybox=True, loc=0,ncol=2)
leg.get_frame().set_alpha(0.5)
plt.xticks(fontsize=14)
plt.xlabel('Epoch', fontsize=15)
plt.ylabel('Objective', fontsize=15)
plt.yticks(fontsize=14)
plt.grid(linestyle='dotted',linewidth=2)

ax = plt.subplot(1, 2, 2)
plt.plot(xaxis, objectiveIEM, label='IEM', marker='^')
plt.plot(xaxis, objectiveEM, label='EM', marker='^')
plt.plot(xaxis, objectiveSAGA, label='FI-EM', marker='^')
plt.plot(xaxis, objectiveoEM_vr, label='SVR-EM', marker='^')
plt.plot(xaxis, objectiveoEM, label='OEM', marker='^')
leg = plt.legend(fontsize=20,fancybox=True, loc=0,ncol=2)
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epoch', fontsize=15)
plt.ylabel('Objective', fontsize=15)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(linestyle='dotted',linewidth=2)
fig.tight_layout()
plt.show()