# Import 3rd Party Software (numpy,ase)
import numpy as np
from ase.calculators.neighborlist import NeighborList

# Function gives the Coordination numbers for all atoms in a cluster given the atoms and generated neighborlist
def cns(neighlist,atoms1):
        cns=[]
        for i in range(len(atoms1)):
		neighs=neighlist.get_neighbors(i)
		neighs=neighs[0].tolist()
		# Because Bothways=True is set on for neighborlists, the identity must be removed from the neighborlists. Idea: alternative would be n-1.
		neighs=[x for x in neighs if x != i]
		cns.append(len(neighs))
	return np.array(cns)
