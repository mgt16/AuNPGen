# Import 3rd party software (ase,numpy)
import numpy as np
from ase.calculators.neighborlist import NeighborList
from ase.calculators.calculator import Calculator, all_changes


class SuttonChen(Calculator):
    implemented_properties = ['energy', 'forces', 'stress']
# Parameters for Gold epsilon (eV), sigma, Angstroms
    default_parameters = {'epsilon': 1.2793e-2,
                          'sigma': 4.33,
                          'rc': 20,
			   'n': 10,
			   'm': 8,
			    'c': 34.408}
    nolabel = True

    def __init__(self, **kwargs):
        Calculator.__init__(self, **kwargs)

    def calculate(self, atoms=None,
                  properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        natoms = len(self.atoms)

        sigma = self.parameters.sigma
        epsilon = self.parameters.epsilon
        rc = self.parameters.rc
	n=self.parameters.n
	m=self.parameters.m
	c=self.parameters.c

        positions = self.atoms.positions
        cell = self.atoms.cell
	
       # Use Neighborlist with bothways=True to calculate Rho values 
        self.nl = NeighborList([rc / 2.] * natoms, self_interaction=False,bothways=True)
	rhos=np.array([])	
        self.nl.update(self.atoms)
	for a1 in range(natoms):
	    neighbors,offsets=self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d**2).sum(1)
            r=r2**0.5
	    ar = sigma/r
            ar[r > rc] = 0.0
	    rho=(ar**(m)).sum()
	    rhos=np.append(rhos,rho)
            
      # uRC is u at cutoff distance 
        uRC = epsilon*(1/2.*(sigma / rc)**n - c*(sigma / rc)**(m/2.))       	
        duRC= epsilon/2.*(n*(sigma/rc**2)**n+c*m/(sigma/rc)**(m/2.)*(sigma/rc**2)**m) 

       # Initialize Outputs
        energy = 0.0
        forces = np.zeros((natoms, 3))
        stress = np.zeros((3, 3))

       # Energy
        for a1 in range(natoms):
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d**2).sum(1)
            r = r2**0.5 
	    ar = sigma/r
            ar[r>rc] = 0.0
	    rho=rhos[a1]
	# Background energy introduced by cutoff Uc(r)=U(r)-U(rc)-dU/dr|rc*(r-rc)
            uR=epsilon*(1/2.*(ar**(n)).sum()-c*(rho)**0.5)
	    uC=uR-uRC-(duRC*(r-rc)).sum() 
            energy += uC

      # Update Neighborlists to not use bothways
        self.nl = NeighborList([rc / 2] * natoms, self_interaction=False)
        self.nl.update(self.atoms)
	# Force and Stress
	for a1 in range(natoms):	
            neighbors, offsets = self.nl.get_neighbors(a1)
            cells = np.dot(offsets, cell)
            d = positions[neighbors] + cells - positions[a1]
            r2 = (d**2).sum(1)
            r = r2**0.5 
	    ar = sigma/r
            ar[r>rc] = 0.0
            fRC=np.array([duRC]*len(r))
	    fRC[r>rc] = 0.0
	# Force definition
	    fks=((-epsilon*(n*(ar)**(n)-c*m/2.*(1/rhos[neighbors]**0.5+1/(rhos[a1]*np.ones(len(neighbors)))**0.5)*ar**(m)))*1./r**2.)[:,np.newaxis]*d
	    f=fks-fRC[:,np.newaxis]*d
	    forces[a1] += f.sum(axis=0) 
            for a2, f2 in zip(neighbors, f):
                forces[a2] -= f2
            stress += np.dot(f.T, d)
	   
        stress += stress.T.copy()
        stress *= -0.5 / self.atoms.get_volume()

        self.results['energy'] = energy
        self.results['forces'] = forces
	self.results['stress'] = stress
