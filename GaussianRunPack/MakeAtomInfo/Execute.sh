#!/bin/bash
python CreateCom.py
bash RunGaussian.sh
python MakeAtomInfo.py > ../AtomInfo.json
