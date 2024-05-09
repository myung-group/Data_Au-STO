import numpy as np

from ase import Atoms
from ase.visualize import view
from ase.io import read, write, Trajectory
from ase.geometry import get_distances
from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii

from scipy.constants import Boltzmann, Planck

from calculate_entropy_enthalpy import *
import subprocess

"""
$ grep THz OUTCAR > frequency.txt
"""

command = "grep THz OUTCAR"
filename = 'frequency.txt'
with open(filename, 'w') as output_file:
    result = subprocess.run(command, shell=True, stdout=output_file)

txt = open(filename)
txt_raw = txt.readlines().copy()
freq_txt = [line.strip().split() for line in txt_raw if not line == '\n']
freq_list = [float(freq_txt[i][-8])*(10**12) for i in range(len(freq_txt))]
v_list = freq_list

T = 298.15
q_v, S_v, H_v = vibrational_contribution(v_list, T)
S = J_per_atom_to_eV_per_atom(S_v+Boltzmann)
TS = J_per_atom_to_eV_per_atom((S_v+Boltzmann)*T)
H = J_per_atom_to_eV_per_atom(H_v+Boltzmann*T)
zpe = zero_point_enery(v_list)
ZPE = J_per_atom_to_eV_per_atom(zpe)

print(f'\n==> ZPE   = {ZPE:9.5f}  eV')
print(f'==> -TS   = {-TS:9.5f}  eV')
print(f'==> H(Cp) = {H:9.5f}  eV\n')
