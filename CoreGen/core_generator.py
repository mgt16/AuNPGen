# Import calculators for Sutton-Chen Potential and Coordination Numbers (cns)
from calculators.sutton_chen import SuttonChen
from calculators.cns import cns

# Import 3rd party software (ase)
from ase.io import read,write
from ase.calculators.neighborlist import NeighborList
from ase.optimize.basin import BasinHopping
from ase.optimize import QuasiNewton
from ase.data import vdw_radii

# Parameters for running basin hopping and coordination numbers
# Number of core structures total to generating using the basin hopping method
nsteps=20
# Atoms object used as the starting structure for basing hopping
au25=read('structures/25_core.xyz')

# Scale is for generating neighborlists and scales vdw_radii to sutton_chen cutoff for defining neighbors
scale=1.1
# Offset defines the acceptable deviation from 'strucutre-dependent stability' model
offset=0.5

# From Model! Establishing weed out parameters for "bad" structures ########
natoms_3=len(au25)**(-1./3)
m=-20.720002254
b=14.23953783
low=m*natoms_3+b-offset
high=m*natoms_3+b+offset
##########################################

# Set SuttonChen calculator and adjust Basin Hopping parameters
au25.set_calculator(SuttonChen())
bh=BasinHopping(atoms=au25,temperature= 10, dr=0.5, optimizer=QuasiNewton)
bh.initialize()
# Generates local_minima.traj file
bh.run(steps=nsteps)

saved=[]
count=0
print low,high
for i in range(nsteps+2):
	# Read in every local minima and write out only ones that match the model
	one=read('local_minima.traj',index='{}'.format(i))
	atomic_numbers=one.get_atomic_numbers()
	# Use vdw_radii and scale to establish cutoffs for defining neighbors
        cutoffs=[scale*vdw_radii[x] for x in atomic_numbers]
	# Define Neighborlist _ bothways=True so that CN calculations easier
	nl=NeighborList(cutoffs,bothways=True)
	nl.update(one)
	coordns=cns(nl,one)
	avgcn=coordns.mean()
	print avgcn
	# Save ones that match the model!
	if avgcn>low and avgcn<high:	
		write('{}.xyz'.format(count),one)
		saved.append(i)
		count+=1

# Convert ones .traj to .xyz for easier viewing
two=read('local_minima.traj',index=':')
write('local_min.xyz',two)
