#!/bin/bash
python CreateCom.py
bash RunGaussian.sh
python MakeNMRInfo.py > ../NMRInfo.json
