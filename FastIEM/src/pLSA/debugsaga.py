
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

#FIRST BATCH EM STEP
oldp = p
EStep()
MStep()
for k in range(0, K):
    for j in range(0, M):
        for i in range(0,N):
            Ssomme += X[i, j] * oldp[i, j, k]
Hsomme = Ssomme
Vsomme = Ssomme
newLoglikelihood = LogLikelihood()
print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " epoch  ", str(newLoglikelihood))
# if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
#     break
objectiveSAGA.append(newLoglikelihood)
oldLoglikelihood = newLoglikelihood

for epoch in range(1, nb_epochs):
    mini_batches_i = [list_indices_i[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
    mini_batches_j = [list_indices_j[epoch][k:k+mini_batch_size] for k in range(0, N, mini_batch_size)]
    for m in range(mini_batch_size):
        # EStep_incremental(mini_batches_i[m])
        # EStep_incremental(mini_batches_j[m])
        EStep()
        for i in mini_batches_i[m]:
            for j in range(0, M):
                denominatorh = 0;
                denominator = 0;
                for k in range(0, K):
                    h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
                    p[i, j, k] = theta[k, j] * lamda[i, k];
                    denominatorh += h[i, j, k];
                    denominator += p[i, j, k];
                if denominatorh == 0:
                    for k in range(0, K):
                        h[i, j, k] = 0;
                else:
                    for k in range(0, K):
                        h[i, j, k] /= denominator;
                        p[i, j, k] /= denominator;
        for k in range(0, K):
            denominator = 0
            for j in range(0, M):
                theta[k, j] = 0
                tmp = 0
                for i in mini_batches_i[m]:
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
                for i in mini_batches_j[m]:
                    h[i, j, k] = listofthetas[i][k, j] * listoflamdas[i, k];
                    tmph +=  (X[i, j] * p[i, j, k] - X[i, j] * h[i, j, k])
                Hsomme[k,j] += tmph
        for i in mini_batches_j[m]:
            listofthetas[i] = theta
            listoflamdas[i,] = lamda[i,]
    newLoglikelihood = LogLikelihood()
    print("[", time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())), "] ", epoch+1, " epoch  ", str(newLoglikelihood))
    # if(oldLoglikelihood != 1 and newLoglikelihood - oldLoglikelihood < threshold):
    #     break
    objectiveSAGA.append(newLoglikelihood)
    oldLoglikelihood = newLoglikelihood