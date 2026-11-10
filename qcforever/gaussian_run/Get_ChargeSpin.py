import sys


def Extract_ChargeSpin(lines):
    """Read the last full Mulliken table in a Gaussian job step."""
    for index in range(len(lines) - 1, -1, -1):
        header = lines[index].strip()
        if header in ("Mulliken charges:", "Mulliken charges and spin densities:"):
            break
    else:
        raise ValueError("Missing Mulliken charges")

    has_spin = header == "Mulliken charges and spin densities:"
    elements, charges, spins = [], [], []
    for line in lines[index + 2:]:
        if line.strip().startswith("Sum of Mulliken charges") and elements:
            return elements, charges, spins
        fields = line.split()
        if len(fields) != 3 + has_spin or fields[0] != str(len(elements) + 1):
            break
        elements.append(fields[1])
        charges.append(float(fields[2]))
        spins.append(float(fields[3]) if has_spin else 0.0)
    raise ValueError("Incomplete Mulliken charge table")


def find_MaxProtic(Atom, Charge):
    HCharge = []
    H_index = []
    for i in range(len(Atom)):
        if Atom[i] == 'H':
            HCharge.append(Charge[i])
            H_index.append(i)
    if H_index != []:
        Max_HCharge = max(HCharge)
        Max_HCharge_index = HCharge.index(Max_HCharge)
        print(Max_HCharge)
        return H_index[Max_HCharge_index]
    else:
        return None
         

if __name__ == '__main__':
    usage = 'Usage; %s jobname' % sys.argv[0]
    try:
        infilename = sys. argv[1]
    except:
        print (usage); sys.exit()
    with open(infilename, 'r') as ifile:
        lines = ifile.readlines()
    Atom, Charge, Spin = Extract_ChargeSpin(lines)
    print(Atom, Charge)
