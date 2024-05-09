import numpy as np

from ase import Atoms
from ase.io import read, write

from ase.constraints import FixAtoms
from ase.calculators.vasp import Vasp

from gofee.candidates import CandidateGenerator, RattleMutation, PermutationMutation, StartGenerator
from gofee.utils import OperationConstraint
from gofee import GOFEE
from gofee.candidates import MakeBox


##1) Set a calculator
calc=Vasp(command='mpirun vasp_gam',xc='beef-vdw',setups='recommended',
          lreal='Auto',algo='A',nelm=200,ediff=1e-5,lwave=False,lcharg=False,
          ismear=0,encut=500,sigma=0.01,ncore=16,luse_vdw=True ,zab_vdw=-1.8867)

##2) Set a molecule which will be positioned to a slab
d = 1.142
molecule = Atoms('CO', positions=[(0, 0, 0), (0, 0, d)])

##3) Set the slab, the molecule will be positioned to this slab
slab = read('POSCAR-template')
c = FixAtoms(indices=np.arange(len(slab)))
slab.pbc=True
slab.set_constraint(c)

##4) Set the box, you can set this parameter using the MakeBox class
stoichiometry = [6]+[8]

box_1 = MakeBox(stoichiometry=stoichiometry,
                slab=slab,
                specified_atoms=79,
                center_point=None,
                bl_factor=2.5,
                shirinkage=None)
box = box_1.make_box()


##5) Set the sg object use the StartGenerator class
"""
If you want to use the method of molecule's random position to the slab, 
you have to set the parameter 'molecule=molecule'

If you want to use the original version, 
you have to set the parameter 'molecule=None' or not set the parameter 'molecule'
"""

radius, spherical_center = box_1.spherical_parameters() # Use the MakeBox class's function

sg = StartGenerator(molecule=molecule,
                    slab=slab, 
                    box=box, 
                    stoichiometry=stoichiometry,
                    radius = radius,
                    center = spherical_center)
#sg.make_mol_structure()

"""
If you want to skip the spherically random generation manner,
just set like below:

sg = StartGenerator(slab, stoichiometry, box, molecule)

do not need to set other parameters, such as course ##5).
"""

##6) Set the box_constraint, the rattle, and the candidate_generator object

# Set up constraint for rattle-mutation
box_constraint = OperationConstraint(box=box)

# initialize rattle-mutation
#n_to_optimize = len(stoichiometry)
n_to_optimize = molecule  # if you want molecule's rattle mutation, you have to set like this
rattle = RattleMutation(n_to_optimize, Nrattle=3, rattle_range=0.3)

candidate_generator = CandidateGenerator([1.0], [sg])



##### Initialize and run search #####
search = GOFEE(calc=calc,
               startgenerator=sg,
               candidate_generator=candidate_generator,
               max_steps=150,
               population_size=10,
               position_constraint=box_constraint)
search.run_new()

