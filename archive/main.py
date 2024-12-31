from ase import Atoms as aseAtoms
from ase import Atom as aseAtom
from ase.io import *
import numpy as np

class Atoms(aseAtoms):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.Assign_radius()
        self.moving_indexes = []
        self.energies = np.array([0.0]*len(self), dtype=np.float64)

    @classmethod
    def read(cls, filename, format=None, symbols=None, **kwargs):
        atoms = read(filename, format=format, **kwargs)
        if symbols:
            unique_symbols = list(set(atoms.symbols))
            if len(unique_symbols) != len(symbols):
                raise ValueError("Number of provided symbols does not match the number of atom types in the file")
            symbol_map = dict(zip(unique_symbols, symbols))
            atoms.symbols = [symbol_map[s] for s in atoms.symbols]
        return cls(atoms)

    def Assign_radius(self, radius):
        # in compliting...
        pass

class 
    

    
    

if __name__ == "__main__":
    target_atoms = Atoms.read("test/min.dump", format="lammps-dump-text", symbols=["B","N"])
    kd_tree = target_atoms.Construct_KDTree()
    target_atoms.assing_radius(1.46)
    
    ion = Atoms()
    ion.append("He", position=(1,0,5))
    ion.velocity = [0,0,-4]
    ion.assing_radius(1.46)





