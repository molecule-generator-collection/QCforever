#Wrapper for Gaussian 16

import glob
import os
import subprocess

from qcforever import laqa_fafoom

hartree2eV = 27.21138602
hartree2kcalmol = 627.503


class g16Object():
    '''Create and handle Gaussian 16 objects.'''
    def __init__(self, sdf_string, gauss_exedir, gauss_scrdir=os.getcwd(),
                 nprocs=1, memory='1GB', jobtype='opt', charge=0, mult=1,
                 qcmethod='pm6', optsteps=200, solvmethod=None, solvent=None,
                 sdf_out='optimized_structures.sdf'):
        """Initialize the g16Object.

        Args(required):
            sdf_string (str): sdf file string
            gauss_exedir (str): environmental variable of ${GAUSS_EXEDIR}
        Args(optional):
            gauss_scrdir (default=`pwd`)
            nprocs       (default=1)
            memory       (default='1GB')
            jobtype      (default='opt')
            charge       (default=0)
            mult         (default=1)
            qcmethod     (default='pm6')
            optsteps     (default=500)
            solvmethod   (default=None)
            solvent      (default=None)
            sdf_out      (default='optimized_structures.sdf')
        Raises:
            KeyError: if the commandline or memory is not defined
        """
        self.sdf_string = sdf_string
        self.gauss_exedir = gauss_exedir
        self.gauss_scrdir = gauss_scrdir
        self.memory = memory
        self.charge = charge
        self.mult = mult
        self.nprocs = nprocs
        self.qcmethod = qcmethod
        self.jobtype = jobtype
        self.optsteps = optsteps
        self.solvmethod = solvmethod
        self.solvent = solvent
        self.sdf_out = sdf_out

    def generate_input(self):
        """Create input files for g16."""
        if self.jobtype == 'gradient':
            jobtype = 'force'
        else:
            jobtype = 'opt=(maxcycle={})'.format(self.optsteps)
        #xyz_string = self.sdf2xyz()
        xyz_string = laqa_fafoom.utilities.sdf2xyz(self.sdf_string)
        coord = xyz_string.split('\n')
        s = '%nprocshared=' + str(self.nprocs) + '\n' + \
            '%mem=' + self.memory + '\n' + \
            '#n ' + self.qcmethod + ' ' + jobtype
        if self.solvmethod is not None:
            s += ' scrf=({}, solvent={})'.format(self.solvmethod, self.solvent)

        s += \
            '\n\ng16 job\n\n' + \
            str(self.charge) + ' ' + str(self.mult) + '\n' + \
            '\n'.join(coord[2:]) + '\n'
        with open('Gau_molecule.com', 'w') as f:
            f.write(s)

    def run_g16(self):
        """Run g16 and write output to 'result.out'. The optimized
        geometry is written to 'self.sdf_string_opt'.

        Warning: this function uses subprocessing to invoke the run.
        The subprocess's shell is set to TRUE.
        Args:
            execution_string (str): e.g. g16 for for parallel version
            /the/complete/path/to/g16
        Raises:
            OSError: if g16_molecule.com not present in the working directory
        """
        success = False
        if os.path.exists('Gau_molecule.com') is False:
            raise OSError("Required input file is not present.")
        os.environ['GAUSS_EXEDIR'] = self.gauss_exedir
        os.environ['GAUSS_SCRDIR'] = self.gauss_scrdir
        print(os.environ['GAUSS_EXEDIR'])
        print(os.environ['GAUSS_SCRDIR'])
        g16 = subprocess.Popen(     \
           "$GAUSS_EXEDIR/g16 Gau_molecule.com", \
            stdout=subprocess.PIPE, shell=True)
        out = subprocess.Popen( \
            ['cat'], stdin=g16.stdout, \
            stdout=open('result.out', 'w'), shell=True) 
        g16.wait()
        out.wait()

        with open('Gau_molecule.log', 'r') as f:
            searchfile = f.readlines()

        opt_conv_key = "Stationary point found"

        not_conv = True
        for line in searchfile:
            if opt_conv_key in line:
                not_conv = False

        #if not_conv:
        #    killfile = open("kill.dat", "w")
        #    killfile.close()

        SCF_conv_key  = "SCF Done"
        coord_bgn_key = "Standard orientation"
        coord_end_key = "Rotational constants (GHZ)"
        grad_bgn_key = "Forces (Hartrees/Bohr)"
        grad_end_key = "Cartesian Forces"

        coord_find = False
        grad_find = False

        for line in searchfile:
            if SCF_conv_key in line:
                list = line.split()
                energy = float(list[4])
                SCF_cycle = int(list[7])
            elif coord_bgn_key in line:
                coord_find = True
                coord = []
            elif coord_end_key in line:
                coord_find = False
            elif coord_find == True:
                list = line.split()
                if len(list) == 6 and list[0] != "Number":
                    coord.append([float(list[3]), float(list[4]), float(list[5])])
            elif grad_bgn_key in line:
                grad_find = True
                grad = []
            elif grad_end_key in line:
                grad_find = False
            elif grad_find == True:
                list = line.split()
                if len(list) == 5 and list[0] != "Number":
                    grad.append([float(list[2]), float(list[3]), float(list[4])])

        self.energy = energy
        self.gradient = grad
        sdf_form = self.sdf_string.splitlines()
        for i in range(4, 4+len(coord)):
            a = sdf_form[i].split()
            c = '{:10.4f}{:10.4f}{:10.4f} {:2}  '.\
                format(coord[i-4][0], coord[i-4][1], coord[i-4][2], a[3])\
                + '  '.join(a[4:])
            sdf_form[i] = c
        self.sdf_string_opt = '\n'.join(sdf_form)
        success = True

        return success

    def get_energy(self, unit='eV'):
        """Get the energy of the molecule.

        Returns:
            energy (float) in eV
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
        """Clean the working direction after the g16 calculation has been
        completed.
        """
        for f in glob.glob("Gau*"):
            os.remove(f)
        for f in glob.glob("fort.*"):
            os.remove(f)
        for f in glob.glob("result.out"):
            os.remove(f)
        for f in glob.glob("core-*"):
            os.remove(f)


def g16_exec(sdf_string, gauss_exedir, gauss_scrdir=os.getcwd(),
             nprocs=1, memory='1GB', jobtype='opt', charge=0, mult=1,
             qcmethod='pm6', optsteps=200, solvmethod=None, solvent=None,
             sdf_out='optimized_structures.sdf'):

    g16_object = g16Object(sdf_string, gauss_exedir, gauss_scrdir, nprocs, memory,
                           jobtype, charge, mult, qcmethod, optsteps,
                           solvmethod, solvent)
    g16_object.clean()
    g16_object.generate_input()
    g16_object.run_g16()
    g16_object.clean()

    energy = g16_object.get_energy()
    grad = g16_object.get_gradient()
    sdf_string_opt = g16_object.get_sdf_string_opt()

    return energy, grad, sdf_string_opt


if __name__ == '__main__':

    sdf_string = '\nH2O\n\n' + \
    '  3  2  0  0  0  0  0  0  0  0999 V2000\n' + \
    '    0.0000    0.0000   -0.3894 O   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '    0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '   -0.7630    0.0000    0.1947 H   0  0  0  0  0  0  0  0  0  0  0  0\n' + \
    '  1  2  1  0  0  0  0\n' + \
    '  1  3  1  0  0  0  0\n' + \
    'M  END\n$$$$'
    gauss_exedir = '/fefs/opt/x86_64/Gaussian/g16'
    gauss_scrdir = os.getcwd()
    nprocs = 8
    memory = '2GB'
    charge = 0
    mult   = 1
    #qcmethod = 'b3lyp/sto-3g'
    qcmethod = 'pm6'
    optsteps = 50
    selvmethod = None
    solvent = None
    solvmethod = 'pcm'
    solvent = 'water'

    print('initial coord')
    print(sdf_string)

    #jobtype = 'gradient'
    jobtype = 'opt'
    energy, grad, sdf_string_opt = g16_exec(sdf_string, gauss_exedir, gauss_scrdir,
                                            nprocs, memory, 
                                            jobtype, charge, mult, qcmethod, optsteps,
                                            solvmethod, solvent)

    print('energy:', energy)
    print('grad')
    for i in range(len(grad)):
        print(grad[i])
    print('opt coord')
    print(sdf_string_opt)
  
