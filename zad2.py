from Bio import PDB
from Bio.PDB.PDBList import PDBList
from Bio.PDB.vectors import Vector
import matplotlib.pyplot as plt

pdbl = PDBList()
fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')
path = input("Path of PDB file:")

pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure("4ywo", path)

def calculate_phi_psi(structure):
    coord_CA = []
    coord_C = []
    coord_N = []
    angles = []
    phi_list = []
    psi_list = []
    for model in structure:
        for chain in model:
            for res in chain:
                for atom in res.get_atoms():
                    if atom.get_name() == 'CA':
                        coord_CA.append(atom.get_coord())
                    if atom.get_name() == 'C':
                        coord_C.append(atom.get_coord())
                    if atom.get_name() == 'N':
                        coord_N.append(atom.get_coord())
    for i in range(0, len(coord_CA)-1):
        if i==0:
            C_prev = Vector(coord_C[len(coord_CA)-1])
        C_prev = Vector(coord_C[i-1])
        N = Vector(coord_N[i])
        CA = Vector(coord_CA[i])
        C = Vector(coord_C[i])
        N_next = Vector(coord_N[i+1])

        phi = PDB.vectors.calc_dihedral(C_prev, N, CA, C)
        psi = PDB.vectors.calc_dihedral(N, CA, C, N_next)
        phi = phi * 180 / 3.14159
        psi = psi * 180 / 3.14159
        phi_list.append(phi)
        psi_list.append(psi)
    angles.extend([phi_list, psi_list])
    return angles

result = calculate_phi_psi(structure)

plt.figure()
plt.scatter(result[0], result[1], s=5)

plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xticks(range(-180, 181, 60))
plt.yticks(range(-180, 181, 60))
plt.axis('scaled')
plt.xlabel("Phi (°)")
plt.ylabel("Psi (°)")

plt.title("Wykres Ramachandrana dla białka 4YWO")

plt.grid(True)
plt.show()