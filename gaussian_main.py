import sys

from qcforever.gaussian_run import GaussianRunPack


def main():
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    #for stargint from gaussian calculation
    #option = "opt pka symm"
    option = "opt homolumo dipole uv symm"
    #option = "symm opt freq nmr uv=/home/sumita/GaussianRun_2.2beta/Samples/UV_peak.dat" 
    #option = "opt homolumo energy dipole deen stable2o2 uv=/home/sumita/GaussianRun_2.2beta/Samples/UV_peak.dat vip vea aip aea" 
    #option = "opt homolumo energy dipole deen stable2o2 uv vip vea aip aea" 
    #option = "opt homolumo energy dipole deen stable2o2 fluor=3" 
    #option = "opt cden homolumo energy dipole deen stable2o2 tadf" 
    #option = "opt nmr uv energy homolumo dipole deen stable2o2 vip vea cden symm"

    test = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 10, option, infilename, solvent='water', restart=False)

    test.para_functional = [0.3]
    test.mem = '5GB'
    test.timexe = 60*60
'''#for geometric constrain
    test.geom_spec = {
        '1 2 3 14': [180.0, 'F'],
        '6 5 4 13': [180.0, 'F'],
        '2 1 7 12': [180.0, 'F'],
        '5 6 8 9' : [180.0, 'F']
    }
'''
'''#Specify spin multiplicity and charge of the target
    test.SpecSpinMulti = 3
    test.SpecTotalCharge = 3
'''
    outdic = test.run_gaussian()
    print (outdic)


if __name__=='__main__':
    main()
