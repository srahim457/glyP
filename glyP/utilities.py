import math
import calc_cp
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

