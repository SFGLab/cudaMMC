#!/usr/bin/python

import os
import argparse
import numpy as np


def convert_to_mmcif(fileName):
    header = """data_3dnome
#
_entry.id cudaMMC
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

    coords = np.loadtxt(fileName, delimiter=' ', skiprows=1)

    i = 1
    for p in coords:
        x, y, z, _ = p

        line = f'ATOM {i} C CA . ALA A 1 {i} ? {x} {y} {z} 1.00 99.99 C\n'
        atoms.append(line)
        i += 1

    with open("%s.cif" % fileName.rstrip(".txt"), 'w') as f:
        f.writelines(atoms)


def main(args):
    indirectory = args.indirectory
    smooth_num = args.smooth_num
    cudaMMC_path = args.cudaMMC_path

    hcms = [os.path.join(indirectory, fname) for fname in os.listdir(indirectory) if fname.endswith('hcm')]
    for hcm in hcms:
        cmd = f"{cudaMMC_path} -a smooth -i {hcm} -r {smooth_num}"
        os.system(cmd)

    txts = [os.path.join(indirectory, fname) for fname in os.listdir(indirectory) if fname.endswith('smooth.txt')]
    for txt in txts:
        convert_to_mmcif(txt)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert to MMCIF using cuMMC")
    parser.add_argument('-i', '--indirectory', required=True, help="Directory containing input files.")
    parser.add_argument('-s', '--smooth_num', required=True,
                        help="Number for smoothing. Number refers to how many nucleotides will be per one bead, i.e. 5000.")
    parser.add_argument('-c', '--cudaMMC_path', required=True,
                        help="Path to the cudaMMC executable.")

    args = parser.parse_args()
    main(args)
