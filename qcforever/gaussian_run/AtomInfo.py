import json
import os

BASE_DIR = os.path.abspath(os.path.dirname(__file__))


def AtomicWeight(Element):
    AtomicWeight = {
        'H': 1.008,'He': 4.003,
        'Li': 6.938, 'Be': 9.012, 'B': 10.81, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'F': 19.00, 'Ne': 20.18,
        'Na': 22.99, 'Mg': 24.80, 'Al':26.98, 'Si':28.08, 'P': 30.97, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.95,
        'K': 39.10, 'Ca': 40.08, 'Sc': 44.96, 'Ti': 47.87, 'V': 50.94, 'Cr': 52.00, 'Mn': 54.94,'Fe': 55.85, 'Co': 58.93,
        'Ni': 58.69, 'Cu': 63.55, 'Zn': 65.38,'Ga': 69.72, 'Ge': 72.63, 'As': 74.92,'Se': 78.97, 'Br': 79.90, 'Kr': 83.80,
        'Rb': 85.47, 'Sr': 87.62, 'Y': 88.91, 'Zr': 91.22, 'Nb': 92.21, 'Mo': 95.95, 'Tc': 99, 'Ru': 101.07, 'Rh': 102.91,
        'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.9, 'Xe': 131.29
        }
    if (Element in AtomicWeight):
        return AtomicWeight[Element]
    else:
        print(f"We don't have the information about {Element}!") 
        exit()


def AtomicNumElec(Element):
    AtomicNumElec = {
        'H' : 1,  'He':  2,
        'Li': 3,  'Be':  4, 'B' : 5,  'C' :  6, 'N' :  7, 'O' :  8, 'F' :  9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P' : 15, 'S' : 16, 'Cl': 17, 'Ar': 18,
        'K':  19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V' : 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27,
        'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y' : 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45,
        'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I' : 53, 'Xe': 54,
        'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63,
        'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 
        'Ta': 73, 'W' : 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81,
        'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 84, 'At': 85, 'Rn': 86
        }
    if isinstance(Element, int):
        if (Element in AtomicNumElec.values()):
            ElementSymbol = [k for k, v in AtomicNumElec.items() if v == Element][0]
            return ElementSymbol
        else:
            print(f"We don't have the information about {Element}!") 
            exit()
    else:  
        if (Element in AtomicNumElec):
            return AtomicNumElec[Element]
        else:
            print(f"We don't have the information about {Element}!") 
            exit()


def ImportDFTdata(filepath, Functional, Basis):
    with open(filepath, 'r') as f:
        data = json.load(f)
    fnc = Functional.upper()
    bas = Basis.upper()
    return data[fnc][bas]


def One_Atom_Energy(Element, Functional, Basis):
    datadir = os.path.join(BASE_DIR, "AtomInfo.json")
    Atom_Energies = ImportDFTdata(datadir, Functional, Basis)
    if (Element in Atom_Energies):
        s = Atom_Energies[Element]
        if s != "none":
            return float(s)
    print(f"We don't have the information about {Element}!") 
    exit()


def One_TMS_refer(Element, Functional, Basis):
    datadir = os.path.join(BASE_DIR, "NMRInfo.json")
    Atom_ref = ImportDFTdata(datadir, Functional, Basis)
    if (Element in Atom_ref):
        s = Atom_ref[Element]
        if s != "none":
            return float(s)
    print(f"We don't have the information about {Element}!") 
    exit()


def O2_MO_refer(Functional, Basis):
    datadir = os.path.join(BASE_DIR, "OxInfo.json")
    O2_MO = ImportDFTdata(datadir, Functional, Basis)
    return float(O2_MO['SOMO_a']), float(O2_MO['LUMO_b'])
