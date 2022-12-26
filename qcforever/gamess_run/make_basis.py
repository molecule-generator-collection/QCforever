import sys

def basis_dissection(basis):

    basis = basis.upper()

    rep_basis = basis.split('-')

    if rep_basis[0] == 'STO':
        GBASIS = 'STO'
        NGAUSS = rep_basis[1].replace('G', '')
        NDFUNC = 0
    else:
        NGAUSS = rep_basis[0]    
        gbasis_dff = rep_basis[1].replace('G', '')
        if '*' in gbasis_dff:
            NDFUNC = gbasis_dff.count('*')
            gbasis = gbasis_dff.replace('*', '')
            GBASIS = 'N'+gbasis
        else:
            GBASIS = 'N'+gbasis_dff
            NDFUNC = 0

    return GBASIS, NGAUSS, NDFUNC

if __name__ == '__main__':
    usage ='Usage; %s infile' % sys.argv[0]

    try:
        in_basis = sys.argv[1]
    except:
        print (usage); sys.exit()

        
    print(basis_dissection(in_basis))


