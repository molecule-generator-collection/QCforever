import sys

from qcforever.gaussian_run import GaussianRunPack


def main():
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        infilename = sys.argv[1]
    except:
        print (usage); sys.exit()

    #for starting from gaussian calculation
    ##For just computing energy
    option = "energy"
    ##For getting the energy after geometry optimization
    #option = "opt energy"
    ##For getting the energy and uv after geometry optimization
    #option = "opt energy uv"
    ##For getting the energy, uv, and fluorescence  after geometry optimization
    #option = "opt energy uv fluor"
    ##For getting the energy, uv, and fluorescence  after geometry optimization with optimization option
    #option = "opt=loose energy uv fluor"
    ##For comparing with uv spectrum data after geometry optimization
    #option = "opt uv=Path2/Samples/UV_peak.dat "
    ##For computing homo/lumo gap and uv after conformation search
    #option = "optconf opt homolumo uv"
    ##For computing homo/lumo gap, dipole moment, decomposition energy, stable to oxygen molecule
    #option = "opt homolumo energy dipole deen stable2o2" 
    ##For computing vertical ionization potential/electronic affinity, adiabatic ionization potential/electronic affinity
    #option = "vip vea aip aea" 
    ##For computing the fluorescence from the third excited state 
    #option = "opt fluor=3" 
    ##For computing the gap between singlet excited state and triplet excited state
    #option = "opt tadf" 

    ##Making the instance with specifying functional, basis set, number of cores.
    #B3LYP/STO-3G with 10 cores
    test = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 10, option, infilename, restart=False)
    #B3LYP/STO-3G with solvent effect after computation chk file will be removed
    #test = GaussianRunPack.GaussianDFTRun('B3LYP', 'STO-3G', 10, option, infilename, solvent='water', restart=False)
    #KTLC-BLYP-BO/3-21G with 10 cores. After computation, the data will be saved as pkl file
    #test = GaussianRunPack.GaussianDFTRun('KTLC-BLYP-BO', '3-21G', 10, option, infilename,  restart=False, pklsave=True)
    #KTLC-BLYP-BO/3-21G with 10 cores. After computation, the data will be saved as pkl file
    #B3LYP/STO-3G after computation chk file will be removed
    #test = GaussianRunPack.GaussianDFTRun('KTLC-BLYP-BO', '3-21G', 20, option, infilename, restart=True, pklsave=True)

    #if you want to set the parameter of LC-DFT.
    #test.para_functional = [0.1961]
    test.mem = '20GB'
    test.timexe = 48*60*60
    #for geometric constrain
    #test.geom_spec = { '1 2 3 14': [180.0, 'F'], '6 5 4 13': [180.0, 'F'], '2 1 7 12': [180.0, 'F'], '5 6 8 9' : [180.0, 'F']}
    #Specify spin multiplicity and charge of the target
    #test.SpecSpinMulti = 3
    #test.SpecTotalCharge = 2
    outdic = test.run_gaussian()

    print (outdic)


if __name__=='__main__':
    main()
