#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
###################
Substructure Search
###################

*Created on Thu Aug 17 15:25 2017 by A. Pahl*

Perform substructure search on the cluster."""

import argparse
import sys
import gzip
import os
import os.path as op
from collections import Counter

from tabulate import tabulate

from rdkit.Chem import AllChem as Chem

# imports for similarity search
# import rdkit.Chem.Scaffolds.MurckoScaffold as MurckoScaffold


def usage():
    print("""Substructure search from multiple Smiles in a tsv file.
    rdkit_sss.py file_to_search file_w_smiles output_dir job_index""")
    sys.exit(1)


def get_value(str_val):
    if not str_val:
        return None
    if isinstance(str_val, list):
        print(str_val)
    try:
        val = float(str_val)
        if "." not in str_val:
            val = int(val)
    except ValueError:
        val = str_val
    return val


def write_rec(f, header, rec):
    line = []
    for h in header:
        line.append(rec.get(h, ""))
    f.write("\t".join(line) + "\n")


def rec_from_cells(header, cells):
    result = {}
    for idx, h in enumerate(header):
        if cells[idx] != "":
            result[h] = cells[idx]
    return result


def iter_csv_mols(fn, header, smiles_col="Canonical_Smiles", max_records=0, tag=True, sep="\t"):
    """Iterate over the mols in a CSV file"""
    if not isinstance(fn, list):
        fn = [fn]
    rec_counter = 0
    for filen in fn:
        if ".gz" in filen:
            f = gzip.open(filen, mode="rt")
        else:
            f = open(filen)
        first_line = True
        for line in f:
            line = line.strip("\n")
            cells = line.split(sep)
            if first_line:
                first_line = False
                header.extend(cells)
                if len(fn) > 1 and tag:
                    header.append("tag")
                continue
            # make a copy with non-empty values
            rec = rec_from_cells(header, cells)  # make a copy with non-empty values
            if smiles_col in rec:
                mol = Chem.MolFromSmiles(rec[smiles_col])
                if not mol:
                    status["Could not generate mol"] += 1
                    continue
            else:
                status["Column {} not found".format(smiles_col)] += 1
            rec["mol"] = mol
            rec_counter += 1
            if max_records > 0 and rec_counter > max_records: break
            if len(fn) > 1 and tag:
                rec["tag"] = filen
            yield rec
        f.close()


def sss_from_file(csv_fn, smiles_fn, output_folder, job_idx, smiles_col="Canonical_Smiles"):
    """Perform a Smiles substructure search.
    The list of Smiles and the CSV file to be searched are loaded from file."""
    smiles_list_orig = open(smiles_fn).read().split("\n")
    smiles_list_orig = [smi.strip() for smi in smiles_list_orig if len(smi) > 2]
    status["Smiles read in"] = len(smiles_list_orig)
    smiles_list = []
    query_mols = []
    for smi in smiles_list_orig:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            smiles_list.append(smi)
            query_mols.append(mol)
        else:
            status["Could not generate query mol"] += 1
    res_fn = op.join(output_folder, "sss_result-{:03}.tsv".format(job_idx))
    result_file = open(res_fn, "w")
    first_line = True
    header_in = []
    hit_ctr = 0
    for mol_ctr, rec in enumerate(iter_csv_mols(csv_fn.format(job_idx), header_in, smiles_col=smiles_col), 1):
        if first_line:
            first_line = False
            header_out = header_in.copy()
            header_out.append("Query")
            result_file.write("\t".join(header_out) + "\n")
            continue
        mol = rec["mol"]
        for idx, qm in enumerate(query_mols):
            if mol.HasSubstructMatch(qm):
                rec["Query"] = smiles_list[idx]
                hit_ctr += 1
                write_rec(result_file, header_out, rec)
                break
    result_file.close()
    status["Molecules searched"] = mol_ctr
    status["Substruct hits"] = hit_ctr


if __name__ == "__main__":
    # file_to_search file_w_smiles output_dir job_index
    parser = argparse.ArgumentParser(
        description="Substructure search from multiple Smiles in a tsv file.")
    parser.add_argument("file_to_search", help="the tsv file to be searched. It needs to contain the formatting place holder for the job id, e.g. '{:03d}'")
    parser.add_argument("file_w_smiles", help="a file with Smiles to be searched, one Smiles per line.")
    parser.add_argument("output_dir", help="where to put the result file.")
    parser.add_argument("job_idx", help="the job index (used for the name of the result file.",
                        type=int)
    parser.add_argument("-s", "--smiles", help="name of the Smiles column in the tsv file.")
    # parser.add_argument("", help="")

    args = parser.parse_args()
    smiles_col = "Canonical_Smiles"
    if args.smiles is not None:
        smiles_col = args.smiles
    err = False
    reason = ""
    if not op.isfile(args.file_to_search.format(args.job_idx)):
        err = True
        reason = "Input file {} could not be found.".format(args.file_to_search)
    elif not op.isfile(args.file_w_smiles):
        err = True
        reason = "File with Smiles {} could not be found.".format(args.file_w_smiles)
    elif not op.isdir(args.output_dir):
        print("Output directory {} does not exist, creating...".format(args.output_dir))
        os.makedirs(args.output_dir, exist_ok=True)

    if err:
        print(reason)
        usage()

    status = Counter()
    sss_from_file(args.file_to_search, args.file_w_smiles, args.output_dir,
                  args.job_idx, smiles_col)

    print(tabulate(status.items()))
    sys.stdout.flush()
