#!/usr/bin/python


import re
import os
import sys
import numpy as np

def convert_to_mmcif(fileName):
    header = """data_3dnome
#
_entry.id 3dnome
#
_audit_conform.dict_name       mmcif_pdbx.dic
_audit_conform.dict_version    5.296
_audit_conform.dict_location   http://mmcif.pdb.org/dictionaries/ascii/mmcif_pdbx.dic
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_asym_id
"""
    atoms = list()
    atoms.append(header)
    with open(fileName, 'r') as f: #open the file
        lines = f.readlines()
        for line in lines:
            values = line.split()
            i = 1
            if(values[4] != "A"):
                i = 0
                num = values[4].split("A")[1]
            else:
                num = values[5]

            new_line = "ATOM " + values[1] + " C " + values[2] + " . " + values[3] + " A  1 " + num + " ? " + values[5+i] + " " + values[6+i] + " " \
            + values[7+i] + " " + values[8+i] + " " + values[9+i] + " A\n"
            atoms.append(new_line)

    with open("%s.mmcif" %fileName.rstrip(".pdb"), 'w') as f:
        f.writelines(atoms)

def convert_to_pdb(fname):
  outf = open("%s.pdb" %fname.rstrip(".txt"), 'w')
  coords = np.loadtxt(fname, delimiter=' ', skiprows=1)

  i = 1
  for p in coords:
    x, y, z, shit = p
    X, Y, Z = round(x,3), round(y,3), round(-z,3)

    if X > 0:
      X = ' '+str(X)+' '
    else:
      X = str(X)+' '
    if len(X) < 8:
      X = X+(9-len(X))*' '
    if Y > 0:
      Y = ' '+str(Y)+' '
    else:
      Y = str(Y)+' '
    if len(Y) < 8:
      Y = Y+(9-len(Y))*' '
    if Z > 0:
      Z = ' '+str(Z)+' '
    else:
      Z = str(Z)+' '
    if len(Z) < 8:
      Z = Z+(9-len(Z))*' '

    I = str(i)
    line = 'ATOM  '+(5-len(I))*' '+I+'   CA ALA A'+(4-len(I))*' '+I+'    '
    line2 = X+Y+Z+'  1.00 99.99\n'
    outf.write(line+line2)
    i = i + 1
  outf.close()


def main():
  indirectory = sys.argv[1]
  smooth_num = sys.argv[2]

  hcms = [os.path.join(indirectory, fname) for fname in os.listdir(indirectory) if fname.endswith('hcm')]
  for hcm in hcms:
    cmd = "./cuMMC -a smooth -i %s -r %s" % (hcm, smooth_num)
    os.system(cmd)

  txts = [os.path.join(indirectory, fname) for fname in os.listdir(indirectory) if fname.endswith('smooth.txt')]
  for txt in txts:
    convert_to_pdb(txt)
    convert_to_mmcif("%s.pdb" %txt.rstrip(".txt"))


if __name__ == '__main__':
  main()
