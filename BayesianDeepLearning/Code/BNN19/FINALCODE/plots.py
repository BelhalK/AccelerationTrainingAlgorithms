import pickle
import matplotlib.pyplot as plt
import numpy as np
import os
print(os.getcwd())




with open ('newlosses/adamloss', 'rb') as fp:
    adam = pickle.load(fp)
nb_epochs = len(adam)
with open ('newlosses/bbbloss', 'rb') as fp:
    bbb = pickle.load(fp)
with open ('newlosses/momentumloss', 'rb') as fp:
    momentum = pickle.load(fp)
with open ('newlosses/dropoutloss', 'rb') as fp:
    dropout = pickle.load(fp)
with open ('newlosses/missoloss', 'rb') as fp:
    misso = pickle.load(fp)
with open ('newlosses/sagloss', 'rb') as fp:
    sag = pickle.load(fp)



nbepochs = 100
iter = len(adam[0:nbepochs])
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16, 4))
axes.set_yscale('log')
#axes.set_xscale('log')
plt.plot(np.arange(iter), adam[0:nbepochs] , label='ADAM', marker='.')
plt.plot(np.arange(iter), bbb[0:nbepochs] , label='BBB', marker='.')
plt.plot(np.arange(iter), momentum[0:nbepochs] , label='Momentum', marker='.')
plt.plot(np.arange(iter), dropout[0:nbepochs] , label='Dropout', marker='.')
plt.plot(np.arange(iter), sag[0:nbepochs] , label='SAG', marker='.')
plt.plot(np.arange(iter), misso[0:nbepochs] , label='MISSO', marker='.')
leg = plt.legend(fontsize=16,fancybox=True, loc=0,ncol=2)
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epochs', fontsize=20)
plt.ylabel('Negated ELBO', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(linestyle='dotted',linewidth=2)
plt.show()

fig.savefig("bnn_with_dropout.png",bbox_inches = 'tight')

