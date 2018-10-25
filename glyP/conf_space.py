import numpy as np
import os
from conformer import *
from copy import copy as cp


class Space(list):

    '''A conformational space consisting of all conformers sfound in specified directory.
    The directory tree should have a structure:
    'molecule'/*/*log
    if directory 'molecule' holds a directory 'experimetnal', an attibute self.expIR is 
    created using the data found there. 
    for different molecules, different lists can (meaning should!) be made.'''

    _temp = 298.15
    _kT=0.0019872036*_temp
    _Ha2kcal=627.5095

    def __init__(self, molecule):

        for (root, dirs, files) in os.walk('./'+molecule):
            for dirname in dirs:
                print dirname
                #oldername = os.path.basename(dirpath)
                if dirname == 'experimental':
                    self.expIR = np.genfromtxt(molecule+'/'+dirname+'/exp.dat')
                for ifiles in os.walk(molecule+'/'+dirname):
                    for filename in ifiles[2]:
                        if filename.endswith('.log'):
                            self.append(Conformer(molecule+'/'+dirname+'/'+filename))

    def __str__(self):
         
        '''Prints a nice table with coded molecular values'''

        print "%20s%20s%20s%20s\n" %('id', 'E', 'H', 'F'),
        for conf in self: 
            print "%20s%20.2f%20.2f%20.2f\n" %(conf._id, conf.E*self._Ha2kcal, conf.H*self._Ha2kcal, conf.F*self._Ha2kcal),
        return ''

    def gaussian_broadening(self, broaden=5):

        ''' Performs gaussian broadening for the set''' 

        for conf in self: conf.gaussian_broadening(broaden)

    def reference_to_zero(self, energy_function='E'):

       '''Finds a conformer with the lowest specified energy function and 
        references remainins conformers to this.'''

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
