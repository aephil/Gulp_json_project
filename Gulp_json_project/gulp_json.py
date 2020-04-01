# coding: utf-8


# In[20]:


import pandas as pd
from netCDF4 import Dataset
import numpy as np
import argparse
import json
import math
from numpy import cos, sin
#from phononweb import Phonon
from .jsonencoder import *
from .gulp_pw import *
pw = gulp_pw()

class gulp_json():
    def __init__(self, filename):
        self.features(filename)
        self.write_json(filename)
        
    def features(self, filename):
        main = open("%s.eig"%filename)
        main_lines = main.readlines()

        #natoms:
        natoms = int(main_lines[0].split()[0])
        self.natoms = natoms

        #"atom_numbers", "car", "atomic_numbers":
        atom_car = []
        for i in range(1, natoms+1):
                atom_car.append((main_lines[i].split()))
        atom_car = np.array(atom_car)
        atom_numbers  = atom_car[:, 0]
        car = atom_car[:, 1:4]
        car = [[float(i) for i in j] for j in car]
        atom_numbers = [int(i) for i in atom_numbers]
        atomic_numbers = np.unique(atom_numbers)
        self.atomic_numbers = atomic_numbers
        self.atom_numbers = atom_numbers
        
        
        #nphons, nqpoints:
        nqpoints = int(main_lines[natoms+1].split()[0])-1
        nphons = int(main_lines[natoms+2].split()[0]) 
        self.nqpoints = nqpoints
        self.nphons = nphons
        
        #reps, chemical formula
        reps = (3, 3, 3)
        self.reps = reps
        
        chemical_formula = pw.get_chemical_formula(atom_numbers)
        self.chemical_formula = chemical_formula

        # In[3]:


        #pos, atom_types, chemical_symbols

        with open(f"{filename}.gin") as inp_txt:
	        inp_lines = inp_txt.readlines()
            
        cell_index = 0
        found = False
        
        for i in inp_lines:
            key = i.split()
            try:
                if key[0][:4].lower() == 'frac':
                    found = True
                    break
            except:
                pass
            cell_index += 1

        if not found:
            raise Exception("No 'fractional' keyword found in GULP input file.")
            
        atom_types = []
        pos = []

        while len(atom_types) < natoms:
            cell_index += 1
            
            try:
                # get next line, without comments
                current_line = inp_lines[cell_index].split('#')[0]
            except IndexError:
                raise Exception(f"Only found {len(atom_types)} atoms (out of {natoms} "
                                "expected) in GULP input file.")
            
            # break into 'words'
            current_line = current_line.split()
            
            # format is type, core/shel, positions * 3, other parameters
            if len(current_line) < 4:
                continue
            try:
                current_pos = [float(f) for f in current_line[1:4]]
            except ValueError:
                if current_line[1] == 'core': # ignore shells for phonons
                    current_pos = [float(f) for f in current_line[2:5]]

            atom_types.append(current_line[0])
            pos.append(current_pos)

        atom_types = np.array(atom_types)
        # not sure why this needs to be an array, but for consistency with existing code

        chemical_symbols = np.unique(atom_types).tolist()
        pos = np.array(pos)
        
        self.pos = pos
        self.atom_types = atom_types
        self.chemical_symbols  = chemical_symbols
    
        
        #cell

        cell_index =0
        found = False
        
        for i in inp_lines:
            key = i.split()
            if key == ['cell']:
                found = True
                break 
            cell_index = cell_index +1

        if not found:
            raise Exception("No 'cell' keyword found in GULP input file. ('vectors' not yet implemented.)")
            
        cell_index = cell_index + 1
        cell_coeff = inp_lines[cell_index].split()
        cell_coeff = [float(i) for i in cell_coeff]
        cell = np.array(pw.lat_vec(cell_coeff))
        self.cell = cell

        # In[4]:
        self.car = pw.red_car(pos, cell)

        #getting the eigenvalues in the right format.

        eig_val_raw = np.loadtxt("%s.disp"%filename)[:, 1]
        eigenvalues = eig_val_raw.reshape([len(eig_val_raw)//nphons, nphons])
        self.eigenvalues = eigenvalues
        
        #getting the qpoints in right format:
        arr = []
        for i in main_lines:
            arr.append(i.split())
        q = []    
        for i in arr:
            if i[0]=='K':
                q.append([float(i[3]), float(i[4]), float(i[5])])
        q = np.array(q)
        qpoints = q[:nqpoints, :]
        self.qpoints = qpoints

        # In[5]:


        #getting the eigenvectors in the correct format


        eig_vec_raw = open("%s.eig"%filename)

        eig_vec_text = open("eig_vec.txt", "w+")

        #eig_vec.seek(0)
        #eig_vec.truncate()
        rd = eig_vec_raw.readlines()[natoms+3: ]
        j = 1
        qpt_temp = nqpoints
        while(qpt_temp):
                for i in range(0, nphons):
                    for l in range(3, 3+natoms):
                        temp = j + i*(2 + natoms) + l
                        eig_vec_text.write(rd[temp-1])
                        #print rd[temp-1]
                j = j + nphons*(natoms + 2) + 1
                qpt_temp = qpt_temp - 1

        eig_vec_text = open("eig_vec.txt")
        eig_vec_line = eig_vec_text.readlines()

        eig_vec1 = open("eig_vec1.txt", "w+")
        for i in range(0, (nphons*natoms)):
            eig_vec1.write(eig_vec_line[i])
        eig_vec1.close()

        eig_vec2 = open("eig_vec2.txt", "w+")
        for i in range(nphons*natoms, (len(eig_vec_line)-(nphons*natoms))):
            eig_vec2.write(eig_vec_line[i])
        eig_vec2.close()


        eig1 = np.loadtxt("eig_vec1.txt")
        zarray = np.zeros([nphons*natoms, 3])
        eig_vector1 = pd.concat([pd.DataFrame(eig1), pd.DataFrame(zarray)], axis =1, ignore_index = True)
        eig_vector2 = pd.DataFrame(np.loadtxt("eig_vec2.txt"))
        eig_vector_full = pd.concat([eig_vector1, eig_vector2], ignore_index = True)
        eig_vector_full = pd.concat([eig_vector_full, eig_vector1], ignore_index = True)
        eig_vector_full.head()


        eig_vec_real = eig_vector_full.iloc[:, 0:3].values
        eig_vec_complex = eig_vector_full.iloc[:, 3:6].values
        eig_vector_mat = np.zeros([nphons*nqpoints*natoms, 3], dtype = complex)
        for i in range(0, nphons*nqpoints*natoms):
            for p in range(0, 3):
                eig_vector_mat[i][p] = complex(eig_vec_real[i][p], eig_vec_complex[i][p])
        eigen_vectors = np.zeros([nqpoints, nphons, natoms, 3], dtype=complex)
        for k in xrange(nqpoints):
            for n in xrange(nphons):
                for i in xrange(natoms):
                    eigen_vectors[k][n][i] = np.array(eig_vector_mat[k*(nphons*natoms)+n*(natoms)+i])
        nqpoints     = len(qpoints)
        nphons       = nphons 
        eigenvectors = eigen_vectors.view(dtype=float).reshape([nqpoints,nphons,natoms, 3 ,2])
        qpoints      = qpoints
        self.nqpoints = nqpoints
        self.nphons = nphons 
        self.eigenvectors = eigenvectors
        self.qpoints = qpoints


        # In[6]:


        #calculate reciprocal lattice
        rec = pw.rec_lat(cell)
        self.rec = rec
        
        #calculate qpoints in the reciprocal lattice
        car_qpoints = pw.red_car(qpoints,rec)
        self.car_qpoints = car_qpoints
        self.distances = []
        distance = 0
        #iterate over qpoints
        for k in range(1,nqpoints):
            self.distances.append(distance);

        #calculate distances
            step = np.linalg.norm(car_qpoints[k]-car_qpoints[k-1])
            distance += step

        #add the last distances
        self.distances.append(distance)





        # In[7]:


        def collinear(a,b,c):
                    """
                    checkkk if three points are collinear
                    """
                    d = [[a[0],a[1],1],
                         [b[0],b[1],1],
                         [c[0],c[1],1]]
                    return np.isclose(np.linalg.det(d),0,atol=1e-5)

        #iterate over qpoints
        self.highsym_qpts = [[0,'']]
        for k in range(1, nqpoints-1):
            #detect high symmetry qpoints
            if not collinear(qpoints[k-1],qpoints[k],qpoints[k+1]):
                self.highsym_qpts.append((k,''))
        #add final k-point
        self.highsym_qpts.append((nqpoints-1,''))



        # In[8]:


        
    def write_json(self, filename):
                """ Write a json file to be read by javascript
                """
                f = open("%s.json"%filename,"w+")


                #create the datastructure to be put on the json file
                data = {"name":         filename,             # name of the material on the website
                        "natoms":       self.natoms,           # number of atoms
                        "lattice":      self.cell,             # lattice vectors (bohr)
                        "atom_types":   self.atom_types,       # atom type for each atom (string)
                        "atom_numbers": self.atom_numbers,     # atom number for each atom (integer)
                        "formula":      self.chemical_formula, # chemical formula
                        "qpoints":      self.qpoints,          # list of point in the reciprocal space
                        "repetitions":  self.reps,             # default value for the repetititions 
                        "atom_pos_car": self.car,   # atomic positions in cartesian coordinates
                        "atom_pos_red": self.pos,              # atomic positions in reduced coordinates
                        "eigenvalues":  self.eigenvalues,      # eigenvalues (in units of cm-1)
                        "distances":    self.distances,        # list distances between the qpoints 
                        "highsym_qpts": self.highsym_qpts,     # list of high symmetry qpoints
                        "vectors":      self.eigenvectors      # eigenvectors
                       }     


                f.write(json.dumps(data,cls=JsonEncoder,indent=1))
                f.close()




