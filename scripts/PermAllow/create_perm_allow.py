#!/usr/bin/python

import sys
import argparse

class atom():
    def __init__(self):
        self.name = ""
        self.index = 0

def is_na(resn):
    if resn[0] in ["R","D"]: # this accounts for all DNA and old style RNA
        return True
    elif len(resn) < 3:      # this accounts for all other nucleotides that are standard (A, A3, AN, A5, etc)
        return True
    else:
        return False
    
def get_atom_exists(name,atomlist):
    for at in atomlist:
        if at.name ==name:
            return True
    return False

def get_atom_idx(name,atomlist):
    for at in atomlist:
        if at.name == name:
            return at.index
    return -1

def get_indices(namelist,atomlist):
    indices = list()
    for name in namelist:
        index = get_atom_idx(name,atomlist)
        indices.append(index)
    return indices

print('''Script to create perm.allow file
(written by K. Roeder)      
================================
Only works for canonical amino and nucleic
acids with AMBER atom names
''')

parser = argparse.ArgumentParser("create_perm_allow.py")
parser.add_argument("input_file", help="PDB file to the parsed")
parser.add_argument("--offset", default=1, type=int, help="Offset if numbering of residues does not start at 1 in pdb file (should be index of first residue)")
args = parser.parse_args()

data = dict()
fname = args.input_file
offset = args.offset

print("Input file: "+ fname)
print("starting index is {}".format(offset))

curr_res = 0

with open(fname, "r") as f:
    for line in f:
        if line[:4] == "ATOM":
            thisatom = atom()
            thisatom.index = int(line[6:11])
            thisatom.name = line[12:16].strip()
            thisresn = line[17:21].strip()
            thisresi = int(line[22:26])
        else:
            continue
        # assume continuous numbering starting at offset
        if thisresi != curr_res + (offset-1):
            curr_res += 1
            data[curr_res] = dict()
            data[curr_res]["name"] = thisresn
            if is_na(thisresn):
                data[curr_res]["type"] = "na"
            else:
                data[curr_res]["type"] = "aa" 
            data[curr_res]["atoms"] = list()
        data[curr_res]["atoms"].append(thisatom)

print("Data for {} residues parsed".format(curr_res))
ntotal = curr_res

groups = list()

for idx in range(1,ntotal+1):
    resn = data[idx]["name"]
    atomlist = data[idx]["atoms"]
    atnames = [at.name for at in atomlist]
    # first take care of nucleic acids
    if (data[idx]["type"]=="na"):
        # first phosphate groups
        if resn[-1:] not in ["N", "5"]:
            groups.append([get_indices(["OP1","OP2"],atomlist),2,0])
        # next C5'
        groups.append([get_indices(["H5'","H5''"],atomlist),2,0])
        # now the nucleobases
        if "A" in resn:
            groups.append([get_indices(["H61","H62"],atomlist),2,0])
        elif "C" in resn:
            groups.append([get_indices(["H41","H42"],atomlist),2,0])
        elif "G" in resn:
            groups.append([get_indices(["H21","H22"],atomlist),2,0])
        elif "T" in resn:
            groups.append([get_indices(["H71","H72","H73"],atomlist),3,0])
        # finally the hydrogens in the sugar at 2'
        if resn[0] == "D":
            groups.append([get_indices(["H2'","H2''"],atomlist),2,0])
    # alternatively, it is an amino acid
    else:
        # deal with capping groups and N termini first
        if resn == "NHE":
            groups.append([get_indices(["HN1","HN2"],atomlist),2,0])       
        # H1 atoms means it's an N terminus, NME or ACE
        elif "H1" in atnames:
            groups.append([get_indices(["H1","H2","H3"],atomlist),3,0])
        # HA2 means it is GLY
        if ("HA2" in atnames) and ("HA3" in atnames):
            groups.append([get_indices(["HA2","HA3"],atomlist),2,0])
        # beta carbon
        if ("HB2" in atnames) and ("HB3" in atnames):
            if ("HB1" in atnames): # this is ala
                groups.append([get_indices(["HB1","HB2","HB3"],atomlist),3,0])
            else:
                groups.append([get_indices(["HB2","HB3"],atomlist),2,0])
        # gamma carbon -> the name of CG is key here, we only deal with the easy residues here
        if ("CG" in atnames):
            if ("HG2" in atnames) and ("HG3" in atnames):
                groups.append([get_indices(["HG2","HG3"],atomlist),2,0])
            #ASP, which is deprotonated
            elif ("OD1" in atnames) and ("OD2" in atnames) and not ("HD2" in atnames):
                groups.append([get_indices(["OD1","OD2"],atomlist),2,0])
        #delta carbon -> same as above, only the simple cases
        if ("CD" in atnames):
            if ("HD2" in atnames) and ("HD3" in atnames):
                groups.append([get_indices(["HD2","HD3"],atomlist),2,0])
            #GLU, which is deprotonated
            elif ("OE1" in atnames) and ("OE2" in atnames) and not ("HE2" in atnames):
                groups.append([get_indices(["OE1","OE2"],atomlist),2,0])
        #epsilon carbon
        if ("CE" in atnames):
            if "MET" in resn:
                groups.append([get_indices(["HE1","HE2","HE3"],atomlist),3,0])
            elif ("HE2" in atnames) and ("HE3" in atnames):
                groups.append([get_indices(["HE2","HE3"],atomlist),2,0])
        #this has dealt with most common groups (mostly CH2 and CH3), now complete with remaining residues
        if "ARG" in resn:
            groups.append([get_indices(["NH1","NH2","HH11","HH21","HH12","HH22"],atomlist),2,2])
            groups.append([get_indices(["HH11","HH12"],atomlist),2,0])
            groups.append([get_indices(["HH21","HH22"],atomlist),2,0])
        if "ASN" in resn:
            groups.append([get_indices(["HD21","HD22"],atomlist),2,0])
        if "GLN" in resn:
            groups.append([get_indices(["HE21","HE22"],atomlist),2,0])
        if "ILE" in resn:
            groups.append([get_indices(["HG21","HG22","HG23"],atomlist),3,0])
            groups.append([get_indices(["HG12","HG13"],atomlist),2,0])
            groups.append([get_indices(["HD11","HD12","HD13"],atomlist),3,0])            
        if "LEU" in resn:
            groups.append([get_indices(["CD1","CD2","HD11","HD23","HD12","HD22","HD13","HD21"],atomlist),2,3])
            groups.append([get_indices(["HD11","HD12","HD13"],atomlist),3,0])            
            groups.append([get_indices(["HD21","HD22","HD23"],atomlist),3,0])            
        if "LYN" in resn:
            groups.append([get_indices(["HZ2","HZ3"],atomlist),2,0])
        if "LYS" in resn:
            groups.append([get_indices(["HZ1","HZ2","HZ3"],atomlist),3,0])
        if "PHE" in resn:
            groups.append([get_indices(["CD1","CD2","CE1","CE2","HD1","HD2","HE1","HE2"],atomlist),2,3])
        if "THR" in resn:
            groups.append([get_indices(["HG21","HG22","HG23"],atomlist),3,0])
        if "TYR" in resn:
            groups.append([get_indices(["CD1","CD2","CE1","CE2","HD1","HD2","HE1","HE2"],atomlist),2,3])
        if "VAL" in resn:
            groups.append([get_indices(["CG1","CG2","HG11","HG23","HG12","HG22","HG13","HG21"],atomlist),2,3])
            groups.append([get_indices(["HG11","HG12","HG13"],atomlist),3,0])            
            groups.append([get_indices(["HG21","HG22","HG23"],atomlist),3,0])
        #finally, we need to add the C termini
        if ("O" in atnames) and ("OXT" in atnames):
            groups.append([get_indices(["O","OXT"],atomlist),2,0])

print("Detected {} permutational groups".format(len(groups)))

for group in groups:
    for idx in group[0]:
        if idx == -1:
            msg = "Error: atom index is -1 after pattern matching: common reasons - wrong offset specified, poorly formatted pdb, or wrong atom names"
            print(msg)
            sys.exit(1)

outf = open("perm.allow","w")

outf.write(str(len(groups))+"\n")
for group in groups:
    indices = group[0]
    nmem = group[1]
    nswap = group[2]
    outf.write(str(nmem)+" "+str(nswap)+"\n")
    string = ""
    for idx in indices:
        string += " " + str(idx)
    outf.write(string+"\n")
outf.close()
print("perm.allow file written")
