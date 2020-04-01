import pandas as pd
import numpy as np
import math
from math import cos, sin


class gulp_pw():  
    
    def __init__(self):
        self.hartree_cm1 = 219474.63
        self.eV          = 27.211396132
        self.Bohr        = 1.88973

        #from ase https://wiki.fysik.dtu.dk/ase/
        self.chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                                'B',  'C',  'N',  'O',  'F',
                                'Ne', 'Na', 'Mg', 'Al', 'Si',
                                'P',  'S',  'Cl', 'Ar', 'K',
                                'Ca', 'Sc', 'Ti', 'V',  'Cr',
                                'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                                'Zn', 'Ga', 'Ge', 'As', 'Se',
                                'Br', 'Kr', 'Rb', 'Sr', 'Y',
                                'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                                'Rh', 'Pd', 'Ag', 'Cd', 'In',
                                'Sn', 'Sb', 'Te', 'I',  'Xe',
                                'Cs', 'Ba', 'La', 'Ce', 'Pr',
                                'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                                'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                                'Yb', 'Lu', 'Hf', 'Ta', 'W',
                                'Re', 'Os', 'Ir', 'Pt', 'Au',
                                'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                                'At', 'Rn', 'Fr', 'Ra', 'Ac',
                                'Th', 'Pa', 'U',  'Np', 'Pu',
                                'Am', 'Cm', 'Bk', 'Cf', 'Es',
                                'Fm', 'Md', 'No', 'Lr']

        self.atomic_numbers = {}
        for Z, symbol in enumerate(self.chemical_symbols):
             self.atomic_numbers[symbol] = Z

    def red_car(self, red,lat):
        """
        Convert reduced coordinates to cartesian
        """
        return np.array([coord[0]*lat[0]+coord[1]*lat[1]+coord[2]*lat[2] for coord in red])

    def car_red(self, car,lat):
        """
        Convert cartesian coordinates to reduced
        """
        return np.array([np.linalg.solve(np.array(lat).T,coord) for coord in car])

    def rec_lat(self,lat):
        """
        Calculate the reciprocal lattice vectors
        """
        a1,a2,a3 = np.array(lat)
        v = np.dot(a1,np.cross(a2,a3))
        b1 = np.cross(a2,a3)/v
        b2 = np.cross(a3,a1)/v
        b3 = np.cross(a1,a2)/v
        return np.array([b1,b2,b3])

    def deg_rad(self, x):
        return (x*(math.pi))/180


    def lat_vec(self, cell_coeff):
        """Calulate the lattice vectors using cell coeff."""
        
        a = cell_coeff[0]
        b = cell_coeff[1]
        c = cell_coeff[2]
        alp = self.deg_rad(cell_coeff[3])
        bet = self.deg_rad(cell_coeff[4])
        gam = self.deg_rad(cell_coeff[5])
        x = [a, 0, 0]
        y = [b*cos(gam), b*sin(gam), 0]
        z = [c*cos(alp)*cos(gam), c*cos(alp)*sin(gam), c*sin(alp)]
        cell = np.array([x, y, z])
        cell = [[round(i, 3) for i in j] for j in cell]
        
        return cell

    def get_chemical_formula(self, atom_numbers):
            """
            from ase https://wiki.fysik.dtu.dk/ase/
            """
            numbers = atom_numbers
            elements = np.unique(numbers)
            symbols = np.array([self.chemical_symbols[e] for e in elements])
            counts = np.array([(numbers == e).sum() for e in elements])

            ind = symbols.argsort()
            symbols = symbols[ind]
            counts = counts[ind]

            if 'H' in symbols:
                i = np.arange(len(symbols))[symbols == 'H']
                symbols = np.insert(np.delete(symbols, i), 0, symbols[i])
                counts = np.insert(np.delete(counts, i), 0, counts[i])
            if 'C' in symbols:
                i = np.arange(len(symbols))[symbols == 'C']
                symbols = np.insert(np.delete(symbols, i), 0, symbols[i])
                counts = np.insert(np.delete(counts, i), 0, counts[i])

            formula = ''
            for s, c in zip(symbols, counts):
                formula += s
                if c > 1:
                    formula += str(c)

            return formula 


    # In[2]:
