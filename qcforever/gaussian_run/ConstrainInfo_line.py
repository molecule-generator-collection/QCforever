import sys


def get_IntCoorddata(intcoord_dict):
    IntCoord_change = []
    IntCoord_value = []
    IntCoord_type = []

    for key in intcoord_dict.keys():
        IntCoord_change.append(key)
        IntCoord_value.append(float(intcoord_dict[key][0]))
        IntCoord_type.append(intcoord_dict[key][1])

    IntCoord_s = ''
    for i in range(len(IntCoord_change)):
        IntCoord_s += IntCoord_change[i] 
        IntCoord_s += f' {IntCoord_value[i]:.1f}' + '\n'
        if IntCoord_type[i] == 'F':
            IntCoord_s += IntCoord_change[i] + ' F\n'
        if IntCoord_type[i] == 'A':
            IntCoord_s += IntCoord_change[i] + ' A\n'

    return IntCoord_s


def read_IntCoorddata(inputdata):
    IntCoord_change = []
    IntCoord_value = []
    IntCoord_type = []

    try:
        with open(inputdata,'r') as infile:
            lines = infile.readlines()
    except:
        print ("No such file!")
        exit()

    for line in lines:
        line_s = line.split()
        if len(line_s) == 4:
            IntCoord_change.append(f'{line_s[0]:4s} {line_s[1]:4s}')
            IntCoord_value.append(float(line_s[2]))
            IntCoord_type.append(line_s[3])
        elif len(line_s) == 5:
            IntCoord_change.append(f'{line_s[0]:4s} {line_s[1]:4s} {line_s[2]:4s}')
            IntCoord_value.append(float(line_s[3]))
            IntCoord_type.append(line_s[4])
        elif len(line_s) == 6:
            IntCoord_change.append(f'{line_s[0]:4s} {line_s[1]:4s} {line_s[2]:4s} {line_s[3]:4s}')
            IntCoord_value.append(float(line_s[4]))
            IntCoord_type.append(line_s[5])
        elif len(line_s) == '':
            pass
        else:
            print ("There are errors in the process of reading the internal coordinate data!")
            exit()

    IntCoord_s = ''
    for i in range(len(IntCoord_change)):
        IntCoord_s += IntCoord_change[i] 
        IntCoord_s += f'{IntCoord_value[i]:.1f} + \n'
        if IntCoord_type[i] == 'F':
            IntCoord_s += IntCoord_change[i] + 'F\n'

    return IntCoord_s


if __name__ == '__main__':

#    usage = 'Usage; %s infilename' % sys.argv[0]
#    try:
#        infilename = sys. argv[1]
#    except:
#        print (usage); sys.exit()

#    print(read_IntCoorddata(infilename))


    intcoord_dict = {
                    '1 2': [8, 'B'], 
                    '1 2 3 14':[180.0, 'F'],
                    '6 5 4 13':[180.0, 'F']
                    }
    print(get_IntCoorddata(intcoord_dict))

