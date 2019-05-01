
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


rho_saga = 0.01
#REINITIALIZE
lamda = np.random.sample([N, K])
theta = np.random.sample([K, M])
p = zeros([N, M, K])
initializeParameters()
oldLoglikelihood = 1
newLoglikelihood = 1

### SAGA EM
objectiveSAGA = []
epoch = 0
Hsomme = zeros([K, M])
Ssomme = zeros([K, M])
EStep()
MStep()
newLoglikelihood = LogLikelihood()
print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " epoch  ", str(newLoglikelihood))
objectiveSAGA.append(newLoglikelihood)
oldLoglikelihood = newLoglikelihood

listofthetas = {}
for i in range(N):
    listofthetas[i] = theta
listoflamdas = lamda

h = p
oldh = p
for k in range(0, K):
    denominator = 0
    for j in range(0, M):
        Hsomme[k, j] = 0
        for i in range(0, N):
            Hsomme[k, j] += X[i, j] * h[i, j, k]

Ssomme = Hsomme
Vsomme = Hsomme



for epoch in range(1, nb_epochs):
    mini_batches_i = [list_indices_i[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
    mini_batches_j = [list_indices_j[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
    for m in range(mini_batch_size):
        SAGAStep(mini_batches_i[m],mini_batches_j[m])
        MStep()
    newLoglikelihood = LogLikelihood() 
    print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " epoch  ", str(newLoglikelihood))
    # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
    #     break
    objectiveSAGA.append(newLoglikelihood)
    oldLoglikelihood = newLoglikelihood


# with open('losses/sagatest', 'wb') as fp: 
#     pickle.dump(objectiveSAGA, fp)
