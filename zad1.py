import os
import math
import matplotlib.pyplot as plt
import numpy as np

from Bio.PDB.PDBList import PDBList

pdbl = PDBList()
if not os.path.exists("hh/pdb2hhb.ent"):
    fetch_pdb = pdbl.retrieve_pdb_file('2HHB', file_format='pdb')
    fetch_pdb = pdbl.retrieve_pdb_file('4ywo', file_format='pdb')

from Bio import PDB

# Create a CIFParser object
mmcif_parser = PDB.MMCIFParser()
pdb_parser = PDB.PDBParser()

structure = pdb_parser.get_structure("2hhb", "hh/pdb2hhb.ent")


# Iterate over structure components
Koordynaty = []
for model in structure:
  for chain in model:
      for res in chain:
        for atom in res.get_atoms():
          if atom.get_name() == 'CA':
            Koordynaty.append(list(atom.get_coord()))
            
print(Koordynaty)
Koordynaty2 = Koordynaty
Macierz = []
for k in Koordynaty:
    tmp = []
    for l in Koordynaty2:
        distance = math.sqrt((k[0]-l[0])**2+(k[1]-l[1])**2+(k[2]-l[2])**2)
        if distance < 8:
            tmp.append(1)
        else:
            tmp.append(0)
    Macierz.append(tmp)
#print(Macierz)

Macierz = np.array(Macierz)

x_coords, y_coords = np.where(Macierz == 1)

# Create the scatterplot
plt.figure(figsize=(8, 8))
plt.scatter(x_coords, y_coords, c='pink', marker='s', label='Odleglosc (<8 Å)')
plt.gca().invert_yaxis()  # Invert the y-axis to match matrix orientation
plt.xlabel("Indeks")
plt.ylabel("Indeks")
plt.title("Mapa kontaktow")
plt.legend()
plt.grid(False)
plt.show()