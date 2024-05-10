import subprocess
from ase import Atoms

from ase.io import read, write
from scipy.constants import Boltzmann, Planck
from entropy_and_enthalpy import *


"""
This is an example about CO2:

1) w0 = 1      # spin multiplicity of the molecule(adsorbate, CO2)
2) T = 298.15  # K 
3) sigma = 2   # rotational symmetry number
4) R : 1.176   # calculated bondlength by geometry relaxation
5) v_list = [70.725*(10**12), 39.643*(10**12), 18.861*(10**12), 18.820*(10**12)] # frequencies list of CO2
6) molecule    # ase.Atoms object of adsorbate
8) P = 101325  # Pa (= 1 atm); the value must be written as Pa
                               Default is 101325 Pa (1 atm)
9) M           # molecular mass; if None, the M will be automatically calculated.
"""

## 1) Set a gas-phase molecule
molecule = read('POSCAR')    # POSCAR file of molecule

## 2) Calculate degrees of freedom
num = 3*len(molecule) - 5    # linear
# num = 3*len(molecule) - 6  # Non-linear

## 3) Get a list of vibrational frequencies
command = "grep THz OUTCAR"
filename = 'frequency.txt'
with open(filename, 'w') as output_file:
    result = subprocess.run(command, shell=True, stdout=output_file)

txt = open(filename)
txt_raw = txt.readlines().copy()
freq_txt = [line.strip().split() for line in txt_raw if not line == '\n']
freq_list = [float(freq_txt[i][-8])*(10**12) for i in range(len(freq_txt))]
v_list = freq_list[:num]

## 4) Calculate total contribution
total_contribution(v_list = v_list,
                   w0 = 1,
                   T = 298.15,
                   molecule = molecule,
                   sigma = 2,
                   P=None)


