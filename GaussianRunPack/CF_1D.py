import numpy as np


def Integral(x, F):
    Int_F = 0
    for i in range(len(F)-1):
        AvFi = (F[i]+F[i+1]) / 2
        Int_F += (x[i+1]-x[i]) * AvFi
    return Int_F


def Corr_func(fx, f, gx = [], g = []):
    """
    Normalization
    Int_CF = Integral(fx, CF)
    CF = CF / Int_CF
    """
    if g == [] and gx == []:
        gx = fx.copy()
        g = f.copy()
    CF = np.zeros(len(f))
    for i in range(len(f)):
        for j in range(len(g)-i):
               CF[i] += f[j] * g[j+i]
    return fx, CF


# This function is not working. Need to check.
def write_TCF(CFx, CF):        

    ofile_CF = open('CF.dat', 'w')

    for i in range(len(CFx)):
        ofile_TC.write('% 10.5f  % 10.5f   % 10.5f \n' % (CFx[i] ,CF[i]))

    ofile_CF.close()
