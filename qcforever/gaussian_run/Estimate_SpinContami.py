import sys


def Estimate_SpinDiff(lines):
    print ("Start evaluating spin contamination")    
    Ideal_SS = 0.0
    Computed_SS = 0.0
    for line in lines:
        if line.find("Charge =") >= 0:
            line_charge_spin = line.split("=")
            given_SpinMulti = float(line_charge_spin[-1])
            TotalS = (given_SpinMulti-1) / 2
            Ideal_SS = TotalS * (TotalS+1)
            continue
        if line.find("<Sx>=") >= 0:
            line_computed_Spin = line.split()
            Computed_SS = float(line_computed_Spin[-3])
            continue
    return Computed_SS, Ideal_SS


if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    print(Estimate_SpinDiff(lines))
