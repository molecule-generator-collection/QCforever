import math
import numpy as np

from scipy.stats import wasserstein_distance
from scipy.signal import correlate
from scipy import integrate


#from qcforever import gaussian_run
from qcforever import util


pi = math.pi


def read_data(inputdata):
    peak = []
    intensity = []
    UV = {}  # not used
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

def get_wasser_vect(ref_x, target_x):
    d = wasserstein_distance(ref_x, target_x)
    return d

def delete_peak(position, intensity):

    A = np.array(position)
    B = np.array(intensity)
    
    A = A[A*B != 0]

    return A 


def smililarity_dissimilarity(ref_UV_peak, ref_UV_int, target_UV_peak, target_UV_int):
    step = 10
    lower = min(min(ref_UV_peak), min(target_UV_peak)) - 5.0
    upper = max(max(ref_UV_peak), max(target_UV_peak)) + 5.0

    ref_x, ref_broaden = broadening(ref_UV_peak, ref_UV_int, lower, upper, step)
    target_x, target_broaden = broadening(target_UV_peak, target_UV_int, lower, upper, step)

    scipyCF_ref = correlate(ref_broaden, ref_broaden,'same')
    scipyCF_target = correlate(target_broaden, target_broaden,'same')
    scipyCrossCF = correlate(ref_broaden, target_broaden,'same')

    Int_ref = integrate.simpson(scipyCF_ref, x=ref_x, dx=step)
    Int_target = integrate.simpson(scipyCF_target, x=ref_x, dx=step)
    Int_Cross = integrate.simpson(scipyCrossCF, x=ref_x, dx=step)

    S = Int_Cross / math.sqrt(Int_ref*Int_target)
    D = Int_ref + Int_target - 2*Int_Cross

    active_peak = delete_peak(target_UV_peak, target_UV_int)
    wass_dist = get_wasser_vect(ref_UV_peak, active_peak)

    print("Wasserstein: ", wass_dist)
    print ("Similaliry: ", S)
    print ("Dissimilaliry: ", D)

    return S, D, wass_dist
