# -*- coding: utf-8 -*-

"""
Global geometry optimization code using
Look Ahead based on Quadratic Approximation (LAQA) method

Coded on Feb 10 2021

@Author: Michio Katouda (RIST) base on the LAQA code
coded by Kei Terayama (Univ. Tokyo, RIKEN AIP, and Kyoto Univ.)
e-mail: katouda@rist.or.jp
"""

import sys
import time
import datetime

from qcforever import laqa_fafoom

def LAQA_confopt_main(param_file, SMILES=""):

    t_laqaopt_bgn = time.time()
    print(f"\nStart LAQA conforation searh job at {datetime.datetime.now()}")

    LAQA_initgeom(param_file, SMILES)
    LAQA_optgeom(param_file)
	
    t_laqaopt_end = time.time()
    print(f"\nFinish LAQA conforation searh job at {datetime.datetime.now()}")
    print(f"Wall time of LAQA conformation optimization job: {t_laqaopt_end - t_laqaopt_bgn:20.2f} sec.")


