import numpy as np

from ase.io import write, read
from ase.calculators.dftb import Dftb
from ase.build import fcc111,fcc100
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.calculators.vasp import Vasp

from gofee.candidates import CandidateGenerator, StartGenerator, RattleMutation, PermutationMutation
from gofee.utils import OperationConstraint
from gofee import GOFEE

### Define calculator ###
calc=Vasp(command='vasp_gpu',xc='PBE',setups='recommended',lreal='Auto',algo='Normal',nelm=200,ediff=1e-5,lwave=False,lcharg=False,ismear=0,ivdw=20,encut=500,sigma=0.01)

### Set up system ###
template = read('./POSCAR-sto110')
write('POSCAR',template)

# Stoichiometry of atoms to be placed
stoichiometry = 50*[79]

## Box for startgenerator and rattle-mutation
k = 0.15  # Shrinkage fraction from each side of the box in v[0] and v[1] directions.
cell = template.get_cell()
# Initialize box with cell
v = np.copy(cell)
# Set height of box
v[2][2] = 15.0
# Shrink box in v[0] and v[1] directions
v[0] *= (1-2*k)
v[1] *= (1-2*k)
# Chose anker point p0 so box in centered in v[0] and v[1] directions.
#z_max_slab = np.max(template.get_positions()[:,2])
z_max_slab = 18.67
p0 = np.array((0, 0, z_max_slab+0.3)) + k*(cell[0]+cell[1])
# Make box
box = [p0, v]

# initialize startgenerator (used to generate initial structures)
sg = StartGenerator(template, stoichiometry, box)

### Set up candidate generation operations ###
# Set up constraint for rattle-mutation
box_constraint = OperationConstraint(box=box)

# initialize rattle-mutation
n_to_optimize = len(stoichiometry)
rattle = RattleMutation(n_to_optimize, Nrattle=3, rattle_range=4)

candidate_generator = CandidateGenerator([0.2, 0.8], [sg, rattle])

### Initialize and run search ###
search = GOFEE(calc=calc,
               startgenerator=sg,
               candidate_generator=candidate_generator,
               max_steps=200,
               population_size=5,
               position_constraint=box_constraint)
search.run(restart='restart.pickl')
