import os, sys, math, re

def Extract_MO(lines):

        AlphaEigenVal = []
        BetaEigenVal = []
#            for line in lines:
        for line in lines:
            if line.find(" basis functions, ") >=0:
                line_StateInfo = line.split()
            #        print (line_StateInfo[0])
                
                NumBasisFunc = int(line_StateInfo[0])

            if line.find(" alpha electrons ") >=0:
                line_StateInfo = line.split()
            #        print (line_StateInfo[0])
            #        print (line_StateInfo[3])
            
                NumAlphaElec = int(line_StateInfo[0])
                NumBetaElec = int(line_StateInfo[3])

            if line.find("Alpha  occ. eigenvalues --") >=0:

                if len(AlphaEigenVal) == NumBasisFunc:
                    AlphaEigenVal=[]
                
                line_removed = line.replace(" Alpha  occ. eigenvalues -- ", "")
#                line_StateInfo = line_removed.split()
                line_StateInfo = re.split('(..........)',line_removed)
                while '' in line_StateInfo:
                    line_StateInfo.remove('')
                if '\n' in line_StateInfo:
                    line_StateInfo.remove('\n')

                for i in range(len(line_StateInfo)):
                    AlphaEigenVal.append(float(line_StateInfo[i]))

            if line.find("Alpha virt. eigenvalues --") >=0:
                line_removed = line.replace(" Alpha virt. eigenvalues -- ", "")
#                line_StateInfo = line_removed.split()
                line_StateInfo = re.split('(..........)',line_removed)
                while '' in line_StateInfo:
                    line_StateInfo.remove('')
                if '\n' in line_StateInfo:
                    line_StateInfo.remove('\n')

                for i in range(len(line_StateInfo)):
                    AlphaEigenVal.append(float(line_StateInfo[i]))

            if line.find("Beta  occ. eigenvalues --") >=0:

                if len(BetaEigenVal) == NumBasisFunc:
                    BetaEigenVal=[]
                
                line_removed = line.replace("  Beta  occ. eigenvalues -- ", "")
#                line_StateInfo = line_removed.split()
                line_StateInfo = re.split('(..........)',line_removed)
                while '' in line_StateInfo:
                    line_StateInfo.remove('')
                if '\n' in line_StateInfo:
                    line_StateInfo.remove('\n')

                for i in range(len(line_StateInfo)):
                    BetaEigenVal.append(float(line_StateInfo[i]))

            if line.find("Beta virt. eigenvalues --") >=0:
                line_removed = line.replace("  Beta virt. eigenvalues -- ", "")
#                line_StateInfo = line_removed.split()
                line_StateInfo = re.split('(..........)',line_removed)
                while '' in line_StateInfo:
                    line_StateInfo.remove('')
                if '\n' in line_StateInfo:
                    line_StateInfo.remove('\n')
                if '\n' in line_StateInfo:
                    line_StateInfo.remove('\n')

                for i in range(len(line_StateInfo)):
                    BetaEigenVal.append(float(line_StateInfo[i]))

        return NumAlphaElec, NumBetaElec, AlphaEigenVal, BetaEigenVal 

if __name__ == '__main__':

    usage = 'Usage; %s jobname' % sys.argv[0]

    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()


    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
       
    print(Extract_MO(lines))




