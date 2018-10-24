import numpy as np
import os
from conformer import *
from copy import copy as cp


class Space(list):


    _temp = 298.15
    _kT=0.0019872036*_temp
    _Ha2kcal=627.5095

    def __init__(self):

        for (dirname, dirs, files) in os.walk('.'):
            for dirname in dirs:
                #oldername = os.path.basename(dirpath)
                if dirname == 'experimental':
                    self.expIR = np.genfromtxt(dirname+'/exp.dat')
                for ifiles in os.walk(dirname):
                    for filename in ifiles[2]:
                        if filename.endswith('.log'):
                            self.append(Conformer(dirname+'/'+filename))

    def __str__(self):
         
        print "%20s%20s%20s%20s\n" %('id', 'E', 'H', 'F'),
        for conf in self: 
            print "%20s%20.2f%20.2f%20.2f\n" %(conf._id, conf.E*self._Ha2kcal, conf.H*self._Ha2kcal, conf.F*self._Ha2kcal),
        return ''

    def gaussian_broadening(self, broaden=5):

        for conf in self: conf.gaussian_broadening(broaden)

    def reference_to_zero(self, energy_function='E'):

        Eref = 0.0 ; Fref = 0.0 ; Href = 0.0 
        for conf in self: 
              if energy_function == 'E' and  conf.E < Eref: 
                    Eref = cp(conf.E) ; Href = cp(conf.H) ; Fref = cp(conf.F)
              elif energy_function == 'H' and  conf.H < Href: 
                    Eref = cp(conf.E) ; Href = cp(conf.H) ; Fref = cp(conf.F)
              elif energy_function == 'F' and  conf.F < Fref: 
                    Eref = cp(conf.E) ; Href = cp(conf.H) ; Fref = cp(conf.F)

        for conf in self: 
              conf.E -= Eref;  conf.H -= Href ;  conf.F -= Fref
