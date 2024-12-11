#Wrapper for xTB

import glob
import os
import subprocess

from qcforever import laqa_fafoom

hartree2eV = 27.21138602
hartree2kcalmol = 627.503


class xTBObject():
    '''Create and handle xTB objects.'''
    def __init__(self, sdf_string, xtb_call, jobtype='opt', 
                 gfn='2', charge=0, mult=1, optsteps=200,
                 solvmethod=None, solvent='water',
                 sdf_out='optimized_structures.sdf'):
                 
        """Initialize the xTBObject.
        Args(required):
            sdf_string (str): sdf file string
            xtb_call   (str): e.g. xTB for for parallel version
                                    /the/complete/path/to/xtb
        Args(optional):
            jobtype    (default='opt')
            gfn        (default=2)
            charge     (default=0)
            mult       (default=1)
            optsteps   (default=500)
            solvmethod (default=None)
            solvent    (default='water')
            sdf_out    (default=optimized_structures.sdf)
        Raises:
            KeyError: if the commandline or memory is not defined
        """
        self.sdf_string = sdf_string
        self.xtb_call = xtb_call
        self.jobtype = jobtype
        self.gfn = gfn
        self.charge = charge
        self.mult = mult
        self.optsteps = optsteps
        self.solvmethod = solvmethod
        self.solvent = solvent
        self.sdf_out = sdf_out

    def generate_input(self):
        """Create input files for xTB."""
#        with open('xtbin.mol', 'w') as f:
#            f.write(self.sdf_string)

        xyz_string = laqa_fafoom.utilities.sdf2xyz(self.sdf_string)
        with open('xtbin.xyz', 'w') as f:
            f.write(xyz_string)

    def run_xtb(self):
        """Run xTB and write output to 'result.out'. The optimized
        geometry is written to 'xtbopt.mol'.

        Warning: this function uses subprocessing to invoke the run.
        The subprocess's shell is set to TRUE.
        Raises:
            OSError: if xtbin.mol not present in the working directory
        """
        #for defining OMP_STACKSIZE
        os.environ["OMP_STACKSIZE"] = "5G"

        success = False
        if os.path.exists('xtbin.xyz') is False:
            raise OSError('Required input file not present.')

        if self.jobtype == 'gradient':
            com_xtb = self.xtb_call \
                      + ' --gfn{:>2} xtbin.xyz --chrg {} --uhf {} --grad'\
                          .format(self.gfn, self.charge, self.mult)
            if self.solvmethod is not None :
                com_xtb += ' --alpb {}'.format(self.solvent)
            xtb = subprocess.Popen(com_xtb, stdout=subprocess.PIPE, shell=True)
            out = subprocess.Popen(['cat'], stdin=xtb.stdout,
                                   stdout=open('result.out', 'w'), shell=True)
            xtb.wait()
            out.wait()

            with open('gradient', 'r') as f:
                searchfile = f.readlines()
                num_atoms = int((len(searchfile) - 3) / 2)
                self.energy = float(searchfile[1].split()[6])
                grad = []
                for i in range(2+num_atoms, 2+2*num_atoms):
                    a = searchfile[i].split()
                    grad.append([float(a[0]), float(a[1]), float(a[2])])
                self.gradient = grad
                success = True

        else:
            com_xtb = self.xtb_call \
                      + ' --gfn{:>2} xtbin.xyz --chrg {} --uhf {} --opt --cycles {}'\
                          .format(self.gfn, self.charge, self.mult, self.optsteps)
            if self.solvmethod is not None :
                com_xtb += ' --alpb {}'.format(self.solvent)
            xtb = subprocess.Popen(com_xtb, stdout=subprocess.PIPE, shell=True)
            out = subprocess.Popen(['cat'], stdin=xtb.stdout,
                                   stdout=open('result.out', 'w'), shell=True)
            xtb.wait()
            out.wait()

            with open('result.out', 'r') as f:
                searchfile = f.readlines()

            opt_conv_key = "GEOMETRY OPTIMIZATION CONVERGED"

            not_conv = True
            for line in searchfile:
                if opt_conv_key in line:
                    not_conv = False

            #if not_conv:
            #    killfile = open("kill.dat", "w")
            #    killfile.close()

            energy_key = "TOTAL ENERGY"

            for line in searchfile:
                if energy_key in line:
                    energy = float(line.split()[3])
            self.energy = energy

            with open('xtbopt.xyz', 'r') as f:
                self.sdf_string_opt = f.read()
                success = True

        return success

    def get_energy(self, unit='eV'):
        """Get the energy of the molecule.

        Returns:
            energy (float) in eV, Hartree, or kcal/mol
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'energy'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            if unit == 'hartree':
                return self.energy
            elif unit == 'kcal':
                return hartree2kcalmol * self.energy
            else: # unit == 'eV':
                return hartree2eV * self.energy

    def get_gradient(self):
        """Get the gradient of the molecule.
        Returns:
            gradient (float) in [Hartree/bohr]
        Raises:
            AttributeError: if energy hasn't been calculated yet
        """
        if not hasattr(self, 'gradient'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.gradient

    def get_sdf_string_opt(self):
        """Get the optimized sdf string.

        Returns:
            optimized sdf string (str)
        Raises:
            AttributeError: if the optimization hasn't been performed yet
        """
        if not hasattr(self, 'sdf_string_opt'):
            raise AttributeError("The calculation wasn't performed yet.")
        else:
            return self.sdf_string_opt

    def save_to_file(self): 
        if os.path.isfile(self.sdf_out):
            f = open(self.sdf_out, "a")
        else:
            f = open(self.sdf_out, "w")
        s = str(self.sdf_string_opt) + '\n' \
            + ">  <Energy>"+'\n' \
            + str(self.energy)+'\n\n' \
            + "$$$$"+'\n'
        f.write(s)
        f.close()

    def clean(self):
        """Clean the working direction after the xtb calculation has been
        completed.
        """
        for f in glob.glob("xtb*"):
            os.remove(f)
        for f in glob.glob("gfn*"):
            os.remove(f)
        for f in glob.glob("wbo"):
            os.remove(f)
        for f in glob.glob("charges"):
            os.remove(f)
        for f in glob.glob("energy"):
            os.remove(f)
        for f in glob.glob("gradient"):
            os.remove(f)
        for f in glob.glob("result.out"):
            os.remove(f)


def xtb_exec(sdf_string, xtb_call, jobtype='opt', gfn='2',
             charge=0, mult=1, optsteps=100, 
             solvmethod=None, solvent='water',
             sdf_out='optimized_structures.sdf'):

    xtb_object = xTBObject(sdf_string, xtb_call, jobtype, 
                           gfn, charge, mult, optsteps,
                           solvmethod, solvent, sdf_out)
    xtb_object.clean()
    xtb_object.generate_input()
    xtb_object.run_xtb()
    if jobtype == 'gradient':
        unit = 'hartree'
        energy = xtb_object.get_energy(unit)
        grad = xtb_object.get_gradient()
        xtb_object.clean() 
        return energy, grad
    else:
        unit = 'hartree'
        energy = xtb_object.get_energy(unit)
        sdf_string_opt = xtb_object.get_sdf_string_opt()
        xtb_object.clean()
        return energy, sdf_string_opt


if __name__ == '__main__':

    sdf_string = '\nH2O\n\n' + \
    '  3  2  0  0  0  0  0  0  0  0999 V2000\n' + \
    '    0.0000    0.0000   -0.3894 O   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '    0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '   -0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '  1  2  1  0  0  0  0\n' + \
    '  1  3  1  0  0  0  0\n' + \
    'M  END\n$$$$'
    xtb_call = '~/anaconda3/bin/xtb'
    gfn = '2' #'1', '2', 'ff'
    charge = 0
    mult = 1
    optsteps = 50
    solvmethod = None # 'alpb'
    solvent = None    # 'water'

    print('initial coord')
    print(sdf_string)

    jobtype='gradient'
    energy, grad = xtb_exec(sdf_string, xtb_call, jobtype,
                            gfn, charge, mult, optsteps, solvmethod, solvent)
                            
    print('energy:', energy)
    print('grad')
    for i in range(len(grad)):
        print(grad[i])

    jobtype='opt'
    energy, xyz_string_opt = xtb_exec(sdf_string, xtb_call, jobtype,
                            gfn, charge, mult, optsteps, solvmethod, solvent)

    print('energy:', energy)
    print('opt coord')
    print(xyz_string_opt)
