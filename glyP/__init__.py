from conformer import *
from conf_space import *
import calc_cp


def _main():

    A154 = Space('Tri_A154')
    A154.gaussian_broadening(broaden=5)
    A154.reference_to_zero(energy_function='F')
    print A154
    #for conf in A154:  conf.plot_ir(plot_exp = True, exp_data = A154.expIR)

if __name__ == '__main__':

    _main()
