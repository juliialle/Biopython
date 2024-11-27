from Bio import PDB
import Bio.PDB
import Bio.PDB.Residue

#                           wczytanie pliku gruboziarnistego
path = input("Working directory:")

pdb_parser = PDB.PDBParser()
structure = pdb_parser.get_structure('430D_cg', path+'/coarse_grain_430D.pdb')

#                       funkcje przygotowujące listy atomów
def make_atoms_list(structure):
    list_coords = []
    side_chain_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                list_atoms = []
                side_chain = []
                for atom in res.get_atoms():
                    if atom.get_name() in ['P', 'C4\'']: # (backbone)
                        list_atoms.append(atom)
                    else:
                        side_chain.append(atom)
                        if res.resname in ['A', 'G']:
                            if atom.get_name()=='N9':
                                list_atoms.append(atom)
                        if res.resname in ['C', 'U']:
                            if atom.get_name()=='N1':
                                list_atoms.append(atom)     
                list_coords.append(list_atoms) # Listę atomów zapisuję do kolejnej listy. Dzięki temu łatwo mogę wybrać atomy z konkretnej reszty, która mnie interesuje.
                side_chain_atoms.append(side_chain)
    return list_coords, side_chain_atoms

def make_template_atoms_list(templates):
    list_coords = []
    for structure in templates:
        for model in structure:
            for chain in model:
                for res in chain:
                    list_atoms = []
                    for atom in res.get_atoms():
                        list_atoms.append(atom)
                    list_coords.append(list_atoms)
    return list_coords

def make_template_residues_list(templates):
    list_resis = []
    for structure in templates:
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.resname in ['A', 'G']:
                        list_atoms = []
                        for atom in res.get_atoms():
                            if atom.get_name() in ['N9', 'C2', 'C6']:
                                list_atoms.append(atom)
                        list_resis.append(list_atoms)
                    elif res.resname in ['C', 'U']:
                        list_atoms = []
                        for atom in res.get_atoms():
                            if atom.get_name() in ['N1', 'C2', 'C4']:
                                list_atoms.append(atom)
                        list_resis.append(list_atoms)                   
    return list_resis

#                           wczytanie szablonów
adenine = pdb_parser.get_structure('ade', path+'/templates/adenine.pdb')
cytosine = pdb_parser.get_structure('cyt', path+'/templates/cytosine.pdb')
guanine = pdb_parser.get_structure('gun', path+'/templates/guanine.pdb')
uracil = pdb_parser.get_structure('ura', path+'/templates/uracil.pdb')
backboneN1 = pdb_parser.get_structure('rib1', path+'/templates/riboseN1.pdb')
backboneN9 = pdb_parser.get_structure('rib9', path+'/templates/riboseN9.pdb')

templates = {
    "A": adenine,
    "C": cytosine,
    "G": guanine,
    "U": uracil,
    "BACKBONE1": backboneN1,
    "BACKBONE9": backboneN9
}

residues_list = [adenine, cytosine, guanine, uracil]

#                           stworzenie nowej struktury
new_structure = Bio.PDB.Structure.Structure('NewStructure')
model = Bio.PDB.Model.Model(0)
new_structure.add(model)
chain = Bio.PDB.Chain.Chain("A")
model.add(chain)

list_backboneN1_grains_template, no_grains_list1 = make_atoms_list(backboneN1)
list_backboneN9_grains_template, no_grains_list2 = make_atoms_list(backboneN9)
list_atoms_template = make_template_atoms_list(templates.values())
list_cg_backbones, list_cg_residues = make_atoms_list(structure)
list_residues_grains_template = make_template_residues_list(residues_list)

#                       superimposing i dodawanie łańcucha
sup = PDB.Superimposer()

res_index = 0

for i in range(len(list_cg_backbones)):
    res = list_cg_backbones[i][0].get_parent() 
    res_index += 1
    atm_index = 0
    
    if res.resname in ['A', 'G']:
        template_backbone = list_atoms_template[5]
        sup.set_atoms(fixed=list_cg_backbones[i], moving=list_backboneN9_grains_template[0])
    elif res.resname in ['C', 'U']:
        template_backbone = list_atoms_template[4]
        sup.set_atoms(fixed=list_cg_backbones[i], moving=list_backboneN1_grains_template[0])
    
    sup.apply(template_backbone)

    residue = Bio.PDB.Residue.Residue((" ", res_index, " "), res.resname, " ")
    chain.add(residue)

    for template_atom in template_backbone:
        if template_atom.get_name() == "N1" or template_atom.get_name() == "N9":
            continue
        atm_index += 1
        atom_name = template_atom.get_name()
        atom_coordinates = template_atom.get_coord()  
        element = atom_name[0]  

    #       tworzenie nowego atomu na podstawie danych z szablonu
        new_atom = Bio.PDB.Atom.Atom(
            name=atom_name,                    # Nazwa atomu
            coord=atom_coordinates,            # Współrzędne atomu (lista [x, y, z])
            bfactor=0.0,                       # Współczynnik temperaturowy (domyślny)
            occupancy=1.0,                     # Zajętość atomu (domyślnie 1.0)
            altloc=" ",                        # Alternatywne lokacje (domyślnie brak)
            fullname=" " + atom_name,          # Pełna nazwa atomu z wyrównaniem
            serial_number=atm_index,            # Unikalny numer atomu
            element=element                    # Symbol pierwiastka
        )
        residue.add(new_atom)

    #       reszty azotowe

    if res.resname == "A":
        template_residue = list_atoms_template[0]
        moving_atoms = list_residues_grains_template[0]
    elif res.resname == "C":
        template_residue = list_atoms_template[1]
        moving_atoms = list_residues_grains_template[1]
    elif res.resname == "G":
        template_residue = list_atoms_template[2]
        moving_atoms = list_residues_grains_template[2]
    elif res.resname == "U":
        template_residue = list_atoms_template[3]
        moving_atoms = list_residues_grains_template[3]  

    sup.set_atoms(fixed=list_cg_residues[i], moving=moving_atoms)
    sup.apply(template_residue)

    for template_atom in template_residue:
        atm_index += 1
        atom_name = template_atom.get_name()
        atom_coordinates = template_atom.get_coord()  
        element = atom_name[0]  

        new_atom = Bio.PDB.Atom.Atom(
            name=atom_name,                    # Nazwa atomu
            coord=atom_coordinates,            # Współrzędne atomu (lista [x, y, z])
            bfactor=0.0,                       # Współczynnik temperaturowy (domyślny)
            occupancy=1.0,                     # Zajętość atomu (domyślnie 1.0)
            altloc=" ",                        # Alternatywne lokacje (domyślnie brak)
            fullname=" " + atom_name,          # Pełna nazwa atomu z wyrównaniem
            serial_number=atm_index,            # Unikalny numer atomu
            element=element                    # Symbol pierwiastka
        )
        residue.add(new_atom)

#               wyświetlenie całej struktury i zapis do pliku
for model in new_structure:
    for chain in model:
        for res in chain:
            for atom in res.get_atoms():
                print(f"Nazwa reszty: {res.get_resname()} | Nazwa atomu: {atom.get_name()} | Koordynaty: {atom.get_coord()}")

io = Bio.PDB.PDBIO()
io.set_structure(new_structure)
with open("new_structure.pdb", "w") as output_file:
    io.save(output_file)

