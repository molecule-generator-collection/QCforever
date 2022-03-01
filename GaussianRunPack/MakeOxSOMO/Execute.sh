#!/bin/bash
python CreateCom.py
./RunGaussian.sh
python MakeOxInfo.py > ./OxInfo.json
