from conformer import *
from conf_space import *
import numpy as np
from utilities import *



def _main():

    A154 = Space('Tri_A154') #Space(); in conf_space.py. A154 is the name of the object that holds a list of all conformers
    A154.gaussian_broadening(broaden=5) # Space class method
    A154.reference_to_zero(energy_function='F') # Space class method
    print A154
    #for conf in A154:  conf.plot_ir(plot_exp = True, exp_data = A154.expIR)

if __name__ == '__main__':

    _main()
