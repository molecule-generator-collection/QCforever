"""Example: python orca_main.py molecule.sdf ["opt energy homolumo dipole"]."""

import argparse

from qcforever.orca_run import OrcaRunPack


def main():
    parser = argparse.ArgumentParser(description='Run QCforever using ORCA 6.1+')
    parser.add_argument('input_file')
    parser.add_argument('options', nargs='?', default='opt energy homolumo dipole')
    parser.add_argument('--functional', default='B3LYP')
    parser.add_argument('--basis', default='def2-SVP')
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--memory', default='')
    parser.add_argument('--solvent', default='0')
    args = parser.parse_args()
    job = OrcaRunPack.OrcaDFTRun(args.functional, args.basis, args.nproc,
                                args.options, args.input_file, solvent=args.solvent, pklsave=True)
    job.mem = args.memory
    result = job.run_orca()
    print(result)
    return 0 if result['log'] == 'normal' else 1


if __name__ == '__main__':
    raise SystemExit(main())
