import numpy as np
from scipy.optimize import minimize
from scipy.stats import bernoulli, binom
from scipy.stats import multivariate_normal as mvn
import random
from itertools import chain
import ipdb


def initialize():
    # initial guesses for parameters
    pis = np.random.random(2)
    pis /= pis.sum()
    mu0 = 3
    mus = np.array([[1,mu0], [-mu0,1]])
    sigmas = np.array([np.eye(2)] * 2)
    return pis, mus, sigmas
    
def em(xs, pis, mus, sigmas, max_iter=100):
    
    estim_mus = []
    n, p = xs.shape
    k = len(pis)

    ll_old = 0
    for iteration in range(max_iter):
        exp_A = []
        exp_B = []
        ll_new = 0

        # E-step
        ws = np.zeros((k, n))
        for j in range(len(mus)):
            for i in range(n):
                ws[j, i] = pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
        ws /= ws.sum(0)

        # M-step
        pis = np.zeros(k)
        for j in range(len(mus)):
            for i in range(n):
                pis[j] += ws[j, i]
        pis /= n

        mus = np.zeros((k, p))
        for j in range(k):
            for i in range(n):
                mus[j] += ws[j, i] * xs[i]
            mus[j] /= ws[j, :].sum()
        #ipdb.set_trace()
        
        estim_mus.append(mus)

        sigmas = np.zeros((k, p, p))
        for j in range(k):
            for i in range(n):
                ys = np.reshape(xs[i]- mus[j], (2,1))
                sigmas[j] += ws[j, i] * np.dot(ys, ys.T)
            sigmas[j] /= ws[j,:].sum()

        # update complete log likelihoood
        ll_new = 0.0
        for i in range(n):
            s = 0
            for j in range(k):
                s += pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
            ll_new += np.log(s)

        ll_old = ll_new

    return ll_new, pis, mus, sigmas,estim_mus

def iem(xs, pis, mus, sigmas, max_iter=100):

    estim_mus = []
    n, p = xs.shape

    index = []
    epochs = int(round(max_iter/n))
    for _ in range(epochs):
        rand = [int(x) for x in np.random.choice(n, n,replace=False)]
        index.append(rand)
    index = list(chain(*index))

    k = len(pis)

    ll_old = 0
    ws = np.zeros((k, n))
    for iteration in range(max_iter):
        exp_A = []
        exp_B = []
        ll_new = 0

        # E-step
        if iteration == 0:
            for j in range(len(mus)):
                for i in range(n):
                    ws[j, i] = pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
        else:
            for j in range(len(mus)):
                ws[j, index[iteration]] = pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[index[iteration]])
        ws /= ws.sum(0)
        # M-step
        
        for j in range(len(mus)):
            for i in range(n):
                pis[j] += ws[j, i]
        pis /= n

        mus = np.zeros((k, p))
        for j in range(k):
            for i in range(n):
                mus[j] += ws[j, i] * xs[i]
            mus[j] /= ws[j, :].sum()
        #ipdb.set_trace()
        
        estim_mus.append(mus)

        sigmas = np.zeros((k, p, p))
        for j in range(k):
            for i in range(n):
                ys = np.reshape(xs[i]- mus[j], (2,1))
                sigmas[j] += ws[j, i] * np.dot(ys, ys.T)
            sigmas[j] /= ws[j,:].sum()

        # update complete log likelihoood
        ll_new = 0.0
        for i in range(n):
            s = 0
            for j in range(k):
                s += pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
            ll_new += np.log(s)

        ll_old = ll_new

    return ll_new, pis, mus, sigmas,estim_mus


def oem(xs, pis, mus, sigmas,stepsize, max_iter=100):

    
    estim_mus = []
    n, p = xs.shape

    index = []
    epochs = int(round(max_iter/n))
    for _ in range(epochs):
        rand = [int(x) for x in np.random.choice(n, n,replace=False)]
        index.append(rand)
    index = list(chain(*index))

    k = len(pis)
    stats_pis = np.zeros(k)
    stats_mus = np.zeros((k, p))
    stats_sigmas = np.zeros((k, p, p))

    ll_old = 0
    ws = np.zeros((k, n))
    for iteration in range(max_iter):
        exp_A = []
        exp_B = []
        ll_new = 0
        # E-step
        if iteration == 0:
            for j in range(len(mus)):
                for i in range(n):
                    ws[j, i] = pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
            ws /= ws.sum(0)
            pis = np.zeros(k)
            for j in range(len(mus)):
                for i in range(n):
                    stats_pis[j] += ws[j, i]
            mus = np.zeros((k, p))
            for j in range(k):
                for i in range(n):
                    stats_mus[j] += ws[j, i] * xs[i]
            sigmas = np.zeros((k, p, p))
            for j in range(k):
                for i in range(n):
                    ys = np.reshape(xs[i]- mus[j], (2,1))
                    stats_sigmas[j] += ws[j, i] * np.dot(ys, ys.T)
        else:
            for j in range(len(mus)):
                new_ws = pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[index[iteration]])
                pis[j] = pis[j] + stepsize[iteration]*(new_ws - pis[j])
                stats_mus[j] = stats_mus[j] + stepsize[iteration]*(ws[j, index[iteration]]* xs[index[iteration]] - stats_mus[j])
                ys = np.reshape(xs[index[iteration]]- mus[j], (2,1))
                stats_sigmas[j] += stats_sigmas[j] + stepsize[iteration]*(ws[j, index[iteration]] * np.dot(ys, ys.T) - stats_sigmas[j])

        
        # M-step
        pis /= n
        for j in range(k):
            mus[j] = stats_mus[j]
            mus[j] /= ws[j, :].sum()
        estim_mus.append(mus)

        for j in range(k):
            sigmas[j] = stats_sigmas[j]
            sigmas[j] /= ws[j,:].sum()

        # update complete log likelihoood
        ll_new = 0.0
        for i in range(n):
            s = 0
            for j in range(k):
                s += pis[j] * mvn(mus[j], sigmas[j]).pdf(xs[i])
            ll_new += np.log(s)
        ll_old = ll_new

    return ll_new, pis, mus, sigmas,estim_mus

