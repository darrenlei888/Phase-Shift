# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 21:31:20 2022

@author: Darren
"""

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt


#Number of particles
N = 20
#Size of the box
L = 40
# density
dens= N/L


#the integers/quantum numbers
def Is(N):
    if N%2 == 0:
        arr1 = np.linspace(1/2, (N-1)/2, int(N/2))
        arr2 = np.negative(np.flip(np.linspace(1/2, (N-1)/2, int(N/2))))
        arr = np.concatenate((arr2, arr1))
    else:
        arr1 = np.linspace(0, (N-1)/2, int((N+1)/2))
        arr2 = np.negative(np.flip(np.linspace(1, (N-1)/2, int(N/2))))
        arr = np.concatenate((arr2, arr1))
    return arr


#Phase shift (inverse square g/r^2)
def θ(x):
    if x==0:
        p=0
    if x>0:
        p = +np.pi*(λ-1)
    if x<0:
        p = -np.pi*(λ-1)
    return p


#Bethe ansatz equations
def BAeqs(x):
    arr =[]
    #for k in c:
    eqs = []
    for i in range(N):
        sums=0
        for j in range(N):
            if j!=i:
                sums += θ(x[i]-x[j])
        eqs.append(x[i]-2*np.pi/L * Is(N)[i] - sums)
    arr.append(eqs)
    return arr

#initial trial
init = 2*np.pi*Is(N)
#the solved k's
#k = optimize.newton_krylov(BAeqs, init, iter=40)
#print(k)

#Total energy
def energy(p):
    return 0.5*np.sum(np.square(p))

#Total energy
def momentum(p):
    return np.sum(p)

#The strength
#c = np.linspace(-2, 2, 10)
#c = 1
#print(energy(k))


#The interacting strength
g=np.linspace(0, 4, 10)
energies = []
#λs = np.linspace(-1/4, 1/4, 50)
for i in g:
    λ= 0.5*(1 + np.sqrt(1 + 4*i))
    k = optimize.newton_krylov(BAeqs, init, iter=100)
    enr = energy(k)
    mom = momentum(k)
    energies.append(enr)
    print(λ,mom/L,enr/L)
    print(k)

plt.figure(dpi=500)
plt.xlabel('The Interacting Strength g')
plt.ylabel('The Total Energy Of the System E')
plt.title('Bethe Ansatz')
plt.plot(g, energies, ms=1, c='r', label= 'N=%d, L =100000'%N)
plt.grid()
plt.legend()
plt.show()
