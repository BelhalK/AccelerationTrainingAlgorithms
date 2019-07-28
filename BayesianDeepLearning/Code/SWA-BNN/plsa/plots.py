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





xaxis = np.arange(nb_epochs)

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(16, 4))
ax = plt.subplot(1, 2, 1)
plt.plot(xaxis, objectiveIEM_1k, label='iEM', marker='^')
plt.plot(xaxis, objectiveEM_1k, label='EM', marker='^')
plt.plot(xaxis, objectiveSAGA_1k, label='fiEM', marker='^')
plt.plot(xaxis, objectiveoEM_vr_1k, label='sEM-VR', marker='^')
plt.plot(xaxis, objectiveoEM_1k, label='sEM', marker='^')
leg = plt.legend(fontsize=14,fancybox=True, loc=0,ncol=2)
leg.get_frame().set_alpha(0.5)
plt.xticks(fontsize=16)
plt.xlabel('Epoch', fontsize=20)
plt.ylabel('Objective', fontsize=20)
plt.yticks(fontsize=16)
plt.grid(linestyle='dotted',linewidth=2)
pylab.ticklabel_format(axis='y',style='sci',scilimits=(1,4))

ax = plt.subplot(1, 2, 2)
plt.plot(xaxis, objectiveIEM, label='iEM', marker='^')
plt.plot(xaxis, objectiveEM, label='EM', marker='^')
plt.plot(xaxis, objectiveSAGA, label='fiEM', marker='^')
plt.plot(xaxis, objectiveoEM_vr, label='sEM-VR', marker='^')
plt.plot(xaxis, objectiveoEM, label='sEM', marker='^')
leg = plt.legend(fontsize=14,fancybox=True, loc=0,ncol=2)
leg.get_frame().set_alpha(0.5)
plt.xlabel('Epoch', fontsize=20)
plt.ylabel('Objective', fontsize=20)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(linestyle='dotted',linewidth=2)
pylab.ticklabel_format(axis='y',style='sci',scilimits=(1,4))
fig.tight_layout()
plt.show()

fig.savefig("plsa2.png",bbox_inches = 'tight')