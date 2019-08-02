import math
import calc_cp
import numpy as np
from scipy import interpolate
from optparse import OptionParser

def get_distance(at1, at2):

    return math.sqrt((at1[0]-at2[0])**2+(at1[1]-at2[1])**2+(at1[2]-at2[2])**2)

def element_symbol(A): 
    periodic_table = { '1' : 'H', '6' : 'C', '7' : 'N', '8' : 'O'}
    return periodic_table[A]

def calculate_ring(xyz, ring_atoms): 

    sorted_atoms = []
    for i in 'O', 'C0', 'C1', 'C2', 'C3', 'C4': sorted_atoms.append(ring_atoms[i])

    phi, psi, R = calc_cp.cp_values(xyz, sorted_atoms) 
    return phi, psi, R

def deriv(spec,h):
   """ calculate first derivative of function 'spec'
       using the central finite difference method up to 6th order,
       and for the first 3 and last 3 grid points the
       forward/backward finite difference method up to 2nd order.
       ...as used in f77-program and suggested by Zanazzi-Jona...
   """ 
   der_spec =[[i[0],0] for i in spec]

   length=len(spec)
   for i in range(3,length-3):
      der_spec[i][1]=(-1*spec[i-3][1]+9*spec[i-2][1]-45*spec[i-1][1]+45*spec[i+1][1]-9*spec[i+2][1]+1*spec[i+3][1])/(60*h)
   for i in range(0,3):
      der_spec[i][1]=(-11*spec[i][1]+18*spec[i+1][1]-9*spec[i+2][1]+2*spec[i+3][1])/(6*h)
   for i in range(length-3,length):
      der_spec[i][1]=(11*spec[i][1]-18*spec[i-1][1]+9*spec[i-2][1]-2*spec[i-3][1])/(6*h)

   return der_spec


def get_range(tspec,espec,w_incr,shift,start,stop):
   """ determine wavenumber range within the comparison between theoretical
       and experimental spectrum is performed (depends on the shift)
   """
   de1=start+shift-espec[0][0]
   if (de1 >= 0 ):
      de1=int((start+shift-espec[0][0])/w_incr+0.00001)
      enstart=de1
      tnstart=int((start-tspec[0][0])/w_incr+0.00001)
   else:
      de1=int((start+shift-espec[0][0])/w_incr-0.00001)
      enstart=0
      tnstart=int((start-tspec[0][0])/w_incr-de1+0.00001)
   de2=stop+shift-espec[-1][0]
   if (de2 <= 0 ):
      de2=int((stop+shift-espec[-1][0])/w_incr-0.00001)
      enstop=len(espec)+de2
      tnstop=len(tspec)+int((stop-tspec[-1][0])/w_incr-0.00001) 
   else:
      de2=int((stop+shift-espec[-1][0])/w_incr+0.00001)
      enstop=len(espec)
      tnstop=len(tspec)+int((stop-tspec[-1][0])/w_incr-de2-0.00001)
   return tnstart, tnstop, enstart, enstop
 

def integrate(integrand,delta):
   """ integrate using the trapezoid method as Zanazzi-Jona suggested and was used in the f77-program...
   """
   integral = 0.5*(integrand[0][1]+integrand[-1][1])   
   for i in range(1,len(integrand)-1):
      integral += integrand[i][1]
   return integral*delta

def ypendry(spec,d1_spec,VI):
   """ calculate the Pendry Y-function: y=l^-1/(l^-2+VI^2) with l=I'/I (logarithmic derivative),
       J.B. Pendry, J. Phys. C: Solid St. Phys. 13 (1980) 937-44
   """
   y=[[i[0],0] for i in spec]

   for i in range(len(spec)):
      if (abs(spec[i][1]) <= 1.E-7):
         if (abs(d1_spec[i][1]) <= 1.E-7):
            y[i][1] = 0 
         else:
            y[i][1] = (spec[i][1]/d1_spec[i][1])/((spec[i][1]/d1_spec[i][1])**2+VI**2)
      else:
         y[i][1] = (d1_spec[i][1]/spec[i][1])/(1+(d1_spec[i][1]/spec[i][1])**2*(VI**2))
   return y














