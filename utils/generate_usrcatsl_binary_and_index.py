"""
Generate USRCAT binary file and accompanying SMILES file for SimilarityLab

By Steven Shave - stevenshave@gmail.com
USRCAT is a 3D orientation and translation invariant 3D molecular similarity
technique.  Originally by Ballester & Richards, I introduced atom typing in
the creation of UFSRAT, which was then published first by Schreyer & Blundell
(Schreyer, A.M., Blundell, T. USRCAT: real-time ultrafast shape recognition
with pharmacophoric constraints. J Cheminform 4, 27 (2012).
https://doi.org/10.1186/1758-2946-4-27). adding a new atom type of 
ring atoms - in-line with the CREDO atom type.
USRCAT is now integrated into RDKit.
This code makes it a bit easier to use in a production environment.
Note, the indexes written to the binary indicating SMILES line positions are
not zero-indexed. 1 is the first line.

The binary writer produces a binary file and SMILES lookup for the binary
entries, as used by the SimilarityLab molecular similarity platform.
"""
import struct
import argparse
import json
import gzip
from rdkit import Chem
from pathlib import Path
from rdkit import RDLogger
import progressbar
from cheminformatics_tools_std_include import *
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT
import io


def usrcat_write_binary(sdf_file_path: Path, gzip_output_binary:bool=False):
    assert sdf_file_path.exists(), "SDF file not found"
    binary_file_path=Path(str(sdf_file_path)+".usrcatsl.bin")
    if gzip_output_binary:
        binary_file_path=Path(str(binary_file_path)+".gz")
    assert not Path(binary_file_path).exists(), "Output binary exists"
    output_binary=open_file_may_be_gzipped(binary_file_path, "wb")
    output_smiles_index=open_file_may_be_gzipped(Path(str(sdf_file_path)+".usrcatsl.smi"), "w")
    bar = progressbar.ProgressBar(
        prefix="Generating binary"
    )
    pos_and_desc_bytes=bytearray(struct.calcsize(usrcat_binary_struct_format_string))
    num_gets=0
    num_good_mols=0
    sdf_reader=None
    gz_compressed_file=None
    if str(sdf_file_path).endswith(".gz"):
        gz_compressed_file=gzip.open(str(sdf_file_path))
        sdf_reader=Chem.ForwardSDMolSupplier(gz_compressed_file)
    else:
        sdf_reader=Chem.SDMolSupplier(str(sdf_file_path))
    
    for mol in sdf_reader:
        num_gets+=1
        if mol is not None:
            if mol.GetNumHeavyAtoms()>2:
                num_good_mols+=1
                usrcat_descriptos=GetUSRCAT(mol)
                struct.pack_into(usrcat_binary_struct_format_string, pos_and_desc_bytes,0, num_good_mols, *usrcat_descriptos)
                # Note we use num_good mols, this means that the first line is #1, not 0 - the smiles lines are not zero-indexed.
                output_binary.write(pos_and_desc_bytes)
                output_smiles_index.write(Chem.MolToSmiles(mol)+" "+mol.GetProp("_Name")+"\n")
            if num_good_mols%1000==0:
                bar.update(num_gets)
    bar.update(num_gets)
    output_binary.close()
    output_smiles_index.close()
    print("Num gets",num_gets)
    print("Num good mols",num_good_mols)

if __name__ == "__main__":
    RDLogger.DisableLog("rdApp.*")
    parser = argparse.ArgumentParser(
        description="USRCAT molecular similarity calculations - binary writer for Similarity Lab, generates USRCAT descriptor file containing line entries to an accompanying smiles lookup file"
    )
    parser.add_argument(
        "query_sdf", help="SDF file"
    )
    parser.add_argument(
        "-compress_output", help="Gzip the resultant binary", action="store_true"
    )
    args = parser.parse_args()

    usrcat_write_binary(Path(args.query_sdf), args.compress_output)
