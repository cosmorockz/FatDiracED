import itertools
import numpy as np
import scipy.special as sp
from scipy import integrate

cutoff = 0
Q = 3
ThetaCut = np.pi

l_max = Q - 0.5 + cutoff
mm = [-l_max + i for i in range(int(2*l_max+1))]

def SLL_normalization(Q, l, m):
    # Normalization constant for the Spherical Landau Level
    # Only depends on the charge Q, Landau level l and angular momentum m
    return np.sqrt(((2*l+1)/(4*np.pi)) * sp.binom(2*l, l-Q)/sp.binom(2*l, l-m))


def SLL(theta, Q, l, m):
    # Spherical Landau Level
    # Q is the charge, l is the Landau level, m is the angular momentum
    # this function outputs the value (complex number) of the SLL wavefunction at a specified angular coordinate
    u = np.cos(theta/2)  # * np.exp(1j*phi/2)
    v = np.sin(theta/2)  # * np.exp(-1j*phi/2)
    pre = SLL_normalization(Q, l, m) * (-1)**(l-m) * \
        v**(Q-m) * u**(Q+m)  # part before the summation
    # part inside the summation
    sum_part = 0
    for s in range(0, int(l-m+1)):
        if (l-Q >= s) and (l+Q >= l-m-s):
            sum_part += (-1)**s * sp.binom(l-Q, s) * sp.binom(l+Q, l-m-s) * \
                (v*v)**(l-Q-s) * (u*u)**s
    wf = pre * sum_part  # total wavefunction
    return wf

def ThreeParticleSixTermsOI(theta,Q,l,m1,m2,m3,m4,m5,m6):
    wf1 = SLL(theta,Q,l,m1)
    wf2 = SLL(theta,Q,l,m2)
    wf3 = SLL(theta,Q,l,m3)
    wf4 = SLL(theta,Q,l,m4)
    wf5 = SLL(theta,Q,l,m5)
    wf6 = SLL(theta,Q,l,m6)
    
    x = 2 * np.pi * np.sin(theta) * (wf1 * wf2 * wf3 * wf4 * wf5 * wf6)
    
    return x

def ThreeParticleInteractionStrength(Q,l,m1,m2,m3,m4,m5,m6,ThetaCut):
    x = integrate.quad(ThreeParticleSixTermsOI,0,ThetaCut,args=(Q,l,m1,m2,m3,m4,m5,m6,))
    
    return x[0]

# All the indices for three-particle interactions at cutoff = 0
ThreIn = []
for m1, m2 in itertools.product(mm,mm):
    for m3,m4 in itertools.product(mm,mm):
        for m5 in mm:
            m6 = m1 + m2 + m3 - m4 - m5
            if abs(m6) <= l_max:
                ThreIn.append((m1,m2,m3,m4,m5,m6))
ThreIn = list(set(ThreIn))

OIlist = []
for am in ThreIn:
    
    m1 = am[0]
    m2 = am[1]
    m3 = am[2]
    m4 = am[3]
    m5 = am[4]
    m6 = am[5]
    OIlist.append(ThreeParticleInteractionStrength(l_max,l_max,m1,m2,m3,m4,m5,m6,ThetaCut))

np.savetxt("ThreeParticleIntegralsQ"+str(Q)+"cutoff"+str(0)+".dat",np.c_[ThreIn,OIlist],fmt="%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.8f")  
    



