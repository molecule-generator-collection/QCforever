import math, os, sys
import copy
#import matplotlib.pyplot as plt
import numpy as np

import GaussianRunPack.CF_1D

pi = math.pi

def read_data(inputdata):

        peak = []
        intensity = []

        UV = {}
        try:
            with open(inputdata,'r') as infile:
                lines = infile.readlines()

            for line in lines:
                line_s = line.split()
                peak.append(float(line_s[0]))
                intensity.append(float(line_s[1]))

            uvdata = [peak, intensity]

        except:
            print ("There are some errors in the process of reading uv data!")
            exit()

        return uvdata


def gaussian(x, sigma, center, intensity):

        y = 0

        for i in range(len(center)):
                y += intensity[i]*math.exp(-math.pow((x-center[i])/sigma,2)/2)/math.sqrt(2*pi*sigma*sigma)

        return y

def lorentian(x, gamma, center, intensity):

        y = 0

        for i in range(len(center)):
                y += intensity[i]/(pi*gamma*(1+math.pow((x-center[i])/gamma,2)))

        return y

def broadening(peak, intensity, lower, upper, step):

        x = np.arange(lower, upper, step)

        y_broaden = [lorentian(j, 25.0, peak, intensity) for j in x]

        return x, y_broaden

def smililarity_dissimilarity(ref_UV_peak, ref_UV_int, target_UV_peak, target_UV_int):

        step = 10

        lower = min(ref_UV_peak)-50.0
        upper = max(ref_UV_peak)+50.0

        ref_x, ref_broaden = broadening(ref_UV_peak, ref_UV_int, lower, upper, step)
        target_x, target_broaden = broadening(target_UV_peak, target_UV_int, lower, upper, step)

        CF_refx, CF_ref = GaussianRunPack.CF_1D.Corr_func(ref_x, ref_broaden)
        CF_targetx, CF_target = GaussianRunPack.CF_1D.Corr_func(target_x, target_broaden)
        CrossCFx, CrossCF = GaussianRunPack.CF_1D.Corr_func(ref_x, ref_broaden, target_x, target_broaden)

        Int_ref = GaussianRunPack.CF_1D.Integral(CF_refx, CF_ref)
        Int_target = GaussianRunPack.CF_1D.Integral(CF_targetx, CF_target)
        Int_Cross = GaussianRunPack.CF_1D.Integral(CrossCFx, CrossCF)

        #plt.plot(ref_x, ref_broaden, label="reference")
        #plt.plot(target_x, target_broaden, label="target")
        #plt.plot(CrossCFx, CrossCF, label="CrossCF")
        #plt.show()

        #print(Int_ref)
        #print(Int_target)
        #print(Int_Cross)

        S = Int_Cross/math.sqrt(Int_ref*Int_target)
        D = Int_ref+Int_target-2*Int_Cross

        print ("Similaliry: ", S)
        print ("Dissimilaliry: ", D)

        return S, D


