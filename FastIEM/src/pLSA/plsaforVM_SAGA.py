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

def MStepsaga(index_i, index_j):
    # update theta
    for i in index_i:
        for j in range(0, M):
            denominator = 0;
            for k in range(0, K):
                h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
                denominator += h[i, j, k];
            if denominator == 0:
                for k in range(0, K):
                    h[i, j, k] = 0;
            else:
                for k in range(0, K):
                    h[i, j, k] /= denominator;
    for k in range(0, K):
        denominator = 0
        for j in range(0, M):
            theta[k, j] = 0
            tmp = 0
            for i in index_i:
                tmp +=  X[i, j] * p[i, j, k] - X[i, j] * h[i, j, k]
            Vsomme = Hsomme[k,j] + N*tmp
            Ssomme[k,j] = (1- rho_saga)*Ssomme[k,j] + rho_saga*Vsomme
            theta[k, j] = Ssomme[k,j]
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
    for j in range(0, M):
        for k in range(0, K):
            tmph = 0
            for i in index_j:
                h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
                tmph +=  (X[i, j] * p[i, j, k] - X[i, j] * h[i, j, k])
            Hsomme[k,j] += tmph
    for i in index_j:
        listofthetas[i] = theta
        listoflamdas[i,] = lamda[i,]


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
    
# set the default params and read the params from cmd
# datasetFilePath = 'dataset2.txt'
datasetFilePath = 'dataset1.txt'
stopwordsFilePath = 'stopwords.dic'
K = 10    # number of topic
nb_epochs = 4
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

mini_batch_size = 4
# lamda[i, j] : p(zj|di)
# lamda = random([N, K])
lamda = np.random.sample([N, K])

# theta[i, j] : p(wj|zi)
# theta = random([K, M])
theta = np.random.sample([K, M])

# p[i, j, k] : p(zk|di,wj)
p = zeros([N, M, K])
initializeParameters()

# EM algorithm
oldLoglikelihood = 1
newLoglikelihood = 1


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

# ## REINITIALIZE
p = zeros([N, M, K])
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

# with open('losses/localiemloss', 'wb') as fp: 
#     pickle.dump(objectiveIEM, fp)

rho_saga = 0.003
## REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])
p = zeros([N, M, K])
initializeParameters()
p = zeros([N, M, K])
h = p
oldLoglikelihood = 1
newLoglikelihood = 1
### SAGA EM
objectiveSAGA = []
Hsomme = zeros([K, M])
Ssomme = zeros([K, M])
Vsomme = 0

listofthetas = {}
for i in range(N):
    listofthetas[i] = theta
listoflamdas = lamda


for epoch in range(0, nb_epochs):
# for epoch in range(0, 2):
    oldp = p
    if epoch == 0:
        EStep()
        MStep()
        for k in range(0, K):
            for j in range(0, M):
                for i in range(0,N):
                    Ssomme += X[i, j] * oldp[i, j, k]
        Hsomme = Ssomme
        Vsomme = Ssomme
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveSAGA.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood
    else:
        mini_batches_i = [list_indices_i[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
        mini_batches_j = [list_indices_j[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
        for m in range(len(mini_batches_i)):
            EStep_incremental(mini_batches_i[m])
            EStep_incremental(mini_batches_j[m])
            MStepsaga(mini_batches_i[m], mini_batches_j[m])
        newLoglikelihood = LogLikelihood()
        print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " iteration  ", str(newLoglikelihood))
        # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
        #     break
        objectiveSAGA.append(newLoglikelihood)
        oldLoglikelihood = newLoglikelihood


if __name__ == '__main__':
    output()
