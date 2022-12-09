import sys


def Extract_MO(lines):
        count = 0
        NumAlphaElec = 0
        NumBetaElec = 0
        AlphaEigenVal = []
        BetaEigenVal = []
        for line in lines:
            if line.find("Number of alpha electrons") >= 0:
                line_alphaelecInfo = line.split()
                NumAlphaElec = int(line_alphaelecInfo[-1])
                continue
            if line.find("Number of beta electrons") >= 0:
                line_betaelecInfo = line.split()
                NumBetaElec = int(line_betaelecInfo[-1])
                continue
            if line.find("Alpha Orbital Energies") >= 0:
                line_alphaorbInfo = line.split()
                Numalphaorb = int(line_alphaorbInfo[-1])
                count += 1
                continue
            if count == 1:
                AlphaEigenVal_info = line.split()
                for i in range(len(AlphaEigenVal_info)):
                    AlphaEigenVal.append(float(AlphaEigenVal_info[i]))
                if len(AlphaEigenVal) == Numalphaorb: 
                    count = 0
            if line.find("Beta Orbital Energies") >= 0:
                line_betaorbInfo = line.split()
                Numbetaorb = int(line_betaorbInfo[-1])
                count += -1
                continue
            if count == -1:
                BetaEignVal_info = line.split()
                for i in range(len(BetaEignVal_info)):
                    BetaEigenVal.append(float(BetaEignVal_info[i]))
                if len(BetaEigenVal) == Numbetaorb:
                    count =0
        return NumAlphaElec, NumBetaElec,  AlphaEigenVal, BetaEigenVal 


if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        infilename = sys.argv[1]
    except:
        print (usage)
        sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print(Extract_MO(lines))
