import os


def functional_para(functional, para):

    iop_line = '\n'

    if functional == 'lc-blyp':
        if len(para) != 1:
            print('Number of parameters is invalid! ignoring parameters except for the first one!\n')
        omega_card = int(para[0]*(10**9))
        str_omega_card = str(omega_card).zfill(10)
        iop_line += 'iop(3/107=' + str_omega_card + ',3/108=' + str_omega_card + ')'

    elif functional =='cam-b3lyp':
        if len(para) != 3:
            print('Number of parameters is invalid! ignoring parameters except for the first one!\n')
        omega_card = int(para[0]*(10**9))
        alpha_card = int(para[1]*(10**9))
        beta_card = int(para[2]*(10**4))
        str_omega_card = str(omega_card).zfill(10)
        str_alpha_card = str(alpha_card).zfill(10)
        str_beta_card = str(beta_card).zfill(5)
        iop_line += 'iop(3/107=' + str_omega_card + ',3/108=' + str_omega_card + ')\n'
        iop_line += 'iop(3/119=' + str_alpha_card + ',3/120=' + str_alpha_card + ')\n'
        iop_line += 'iop(3/130=' + str_beta_card + ',3/131=' + str_beta_card + ')'

    else:
        print ('Invalid combination between functional and parameters!')

    return iop_line

if __name__ == '__main__':
    print(functional_para('lc-blyp', [0.33]))
    print(functional_para('cam-b3lyp', [0.33,0.33, 0.25]))

