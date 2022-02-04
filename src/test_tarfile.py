import tarfile

# fname = "../data/test_tarring.tar.gz"
# tar = tarfile.open(fname, "r:gz")
# names = tar.getnames()
# print(names)
# file = tar.extractfile("PDBs/pdb3hkt.ent")
# print("Contents: ")
# print(file.read())

# tar.extractall()
# tar.close()

PDB_ID_lower = '4PTH'
PDB_file = open(f"/Volumes/Seagate_Extension_Plus/PDBs/pdb{PDB_ID_lower}.ent", "r")
print(type(PDB_file))