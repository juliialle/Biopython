from Bio.PDB import PDBParser, PDBIO, Select

# wybor atomow z reprezentacją gruboziarnistą
class CoarseGrainSelect(Select):
    def __init__(self, start_processing=False):
        super().__init__()
        # atomy gruboziarniste
        self.backbone_atoms = {"P", "C4'"}
        self.purine_atoms = {"N9", "C2", "C6"}
        self.pyrimidine_atoms = {"N1", "C2", "C4"}
        self.start_processing = start_processing  # flaga zeby rozpocząć przetwarzanie od pierwszego atomu P

    def accept_atom(self, atom):
        atom_name = atom.get_name()

        # tylko ATOM
        if atom.get_parent().get_id()[0].strip() != "":
            return False

        # rozpoczynanie przetwarzania, gdy znajdziemy pierwszy atom P
        if not self.start_processing:
            if atom_name == "P":
                self.start_processing = True
            else:
                return False  # ignoruj wszystko do pierwszego P

        # czy atom należy do backbone, puryn, czy pirymidyn
        if atom_name in self.backbone_atoms:
            return True
        elif atom.get_parent().get_resname() in {"A", "G"} and atom_name in self.purine_atoms:  # Puryny
            return True
        elif atom.get_parent().get_resname() in {"C", "U", "T"} and atom_name in self.pyrimidine_atoms:  # Pirymidyny
            return True
        return False

# konwersja pełnoatomowej struktury na reprezentację gruboziarnistą
def convert_to_coarse_grain(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("RNA", input_pdb)

    # zapisanie struktury i filtrowanie atomow
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=CoarseGrainSelect())

if __name__ == "__main__":
    input_pdb = "430D.pdb"
    output_pdb = "coarse_grain_430D.pdb"

    convert_to_coarse_grain(input_pdb, output_pdb)
    print(f"Gruboziarnista reprezentacja zapisana w pliku: {output_pdb}")
