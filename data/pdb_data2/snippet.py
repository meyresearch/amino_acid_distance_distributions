
import glob
from MDAnalysis import *

names_xray = glob.glob("/xray_data/*.pdb")
names_nmr = glob.glob("nmr_data/*.pdb")
for i in range(len(names_xray)):
    pdb_file = names[i]
    u = Universe(pdb_file)
    calphas = u.select_atoms("name CA and segid "+u.residues.segments.segids[0])
    #Then do the stuff you did before

for i in range(len(names_nmr)):
    pdb_file = names[i]
    u = Universe(pdb_file)
    calphas = u.select_atoms("name CA and segid "+u.residues.segments.segids[0])
    #Then do the stuff you did before