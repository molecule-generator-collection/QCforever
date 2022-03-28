import sys, math
import copy
from numpy import *

def Integral(x, F):

    Int_F = 0

    for i in range(len(F)-1):
        AvFi = (F[i] + F[i+1])/2
        Int_F += (x[i+1]-x[i])*AvFi

    return Int_F

def Corr_func(fx, f, gx = [], g = []):
    
    if g == [] and gx == []:
        gx = fx.copy()
        g = f.copy()

    CF = zeros(len(f))

    for i in range(len(f)):
        for j in range(len(g)-i):
               CF[i] += f[j]*g[j+i]

####Normailization##############
#    Int_CF = Integral(fx, CF)
#
#    CF = CF/Int_CF
################################

    return fx, CF

def write_TCF(CFx, CF):        

    ofile_CF = open('CF.dat', 'w')

    for i in range(len(CFx)):
        ofile_TC.write('% 10.5f  % 10.5f   % 10.5f \n' % (CFx[i] ,CF[i]))

    ofile_CF.close()
