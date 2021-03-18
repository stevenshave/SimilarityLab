"""
By Steven Shave - stevenshave@gmail.com
"""

from base64 import b64encode
import gzip
import io
from pathlib import Path
from rdkit import Chem
import numpy as np
from rdkit.Chem.MolStandardize import rdMolStandardize
import struct

# The usrcat binary format:
# Starts with an unsigned long (L) denoting the molecule start position in the input SDF
# 5 sets of 12 floats (f), denoting 5 USR distributions (all atoms,
# hydrophobic, HBA, HBA, and ring distributions) comprising mean, 
# variance and skew for measurements to P1, P2, P3 and P4.
# Note, the unsigned long for the very first molecule denotes the number of
# molecules in the binary descriptor file.  This is useful when dealing with a
# compressed descriptor file, and seeking to the end would require
# decompression of everything.
usrcat_binary_struct_format_string = "L" + ("f" * 12) * 5

allowed_atom_types = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "B", "P", "Si", "H"]
lfc = rdMolStandardize.LargestFragmentChooser()


def standardise_smiles_remove_salts(smiles):
    smiles = Chem.MolToSmiles(rdMolStandardize.Cleanup(Chem.MolFromSmiles(smiles)))
    if "." in smiles:
        mol = lfc.choose(Chem.MolFromSmiles(smiles))
        smiles = Chem.MolToSmiles(rdMolStandardize.Cleanup(mol))
    return smiles


def standardise_mol_remove_salts(molecule):
    cleaned_molecule = rdMolStandardize.Cleanup(molecule)
    smiles = Chem.MolToSmiles(cleaned_molecule)
    if "." in smiles:
        cleaned_molecule = rdMolStandardize.Cleanup(lfc.choose(cleaned_molecule))
    return cleaned_molecule


def standardise_mol_remove_salts_return_smiles(molecule):
    cleaned_molecule = rdMolStandardize.Cleanup(molecule)
    smiles = Chem.MolToSmiles(cleaned_molecule)
    if "." in smiles:
        smiles = Chem.MolToSmiles(
            rdMolStandardize.Cleanup(lfc.choose(cleaned_molecule))
        )
    return smiles


def open_file_may_be_gzipped(filepath: Path, mode: str = "rt", buffer_size=None):
    """Open a file, if it is gzipped, then open with gzip.

    Open a file and use gzip if necessary.  A word of caution if
    opening an SDF for RDKit's SDMolSupplier, it must be opened
    in binary mode even if it is text

    Args:
        file (Path): File to open as a pathlib.Path
        mode (str, optional): Mode the file should be opened in. Defaults to "rt":str.

    Returns:
        [File object]: Open file object
    """
    if buffer_size is None:
        if filepath.suffix == ".gz":
            return gzip.open(filepath, mode)
        return open(filepath, mode)
    else:
        if filepath.suffix == ".gz":
            return gzip.open(filepath, mode, buffering=buffer_size)
        return open(filepath, mode, buffering=buffer_size)


class SDFReader_WithFilePos:
    """
    Class to read an SDF, returning molecules and their position in the file

    Unfortunately, using RDKit's SDMolSupplier does not return accurate file
    positions when running tell() on the open file.  This is because iterating
    over lines (SDMolSupplier does this) in Python uses a hidden read-ahead
    buffer for improved performance which means that tell() does not return
    the correct file position for the end of reads.  This object solves this
    at the cost of performance by reading the SDF line by line (readline calls
    do not suffer from this hidden buffer, whilst readlines does).

    Use the class like this:
    while sdf_reader.started_reading_at!=sdf_reader.end_position:
        mol, mol_position=sdf_reader.get()
        if mol is not None:
            # DO WORK HERE
    """

    file = None
    end_position = None
    started_reading_at = None
    num_mols_read = 0

    def __init__(self, sdf_filename: Path):
        assert sdf_filename.exists(), "Could not open SDF file."
        # We open the file in binary mode as tell on Windows with unix style
        # line endings can give strange results.
        self.file = open_file_may_be_gzipped(sdf_filename, "rb")
        self.file.seek(0, io.SEEK_END)
        self.end_position = self.file.tell()
        self.file.seek(0)
        self.started_reading_at = self.file.tell()
        self.line = ""

    def get_end_file_pos(self):
        return self.end_position

    def get(self):
        if self.file is None:
            return None, None

        line = self.file.readline().decode("ascii")
        molecule_text = ""
        molecule_text += line
        while (
            line is not ""
            and self.file is not None
            and self.file.tell() != self.started_reading_at
        ):
            line = self.file.readline().decode("ascii")
            if not line:
                return None, None
            molecule_text += line
            if line.startswith("$$$$"):
                self.num_mols_read += 1
                current_mol_file_offset = self.started_reading_at
                self.started_reading_at = self.file.tell()
                return Chem.MolFromMolBlock(molecule_text), current_mol_file_offset
        return None, None

    def get_from_filepos(self, filepos):
        self.file.seek(filepos)
        return self.get()[0]


class KeepNHighest:
    """
    Keep the best N items

    Object allows insertion of object and score, keeping only the best top N scores.
    """

    num_to_keep = None
    cutoff = None
    scores = []
    items = []  # Items contains tuples of (score, ob)

    def __init__(self, num_to_keep):
        self.num_to_keep = num_to_keep
        self.cutoff = float("-inf")

    def insert(self, item, score):
        if score > self.cutoff:
            insertion_pos = np.searchsorted(-np.array(self.scores), -score)
            self.scores.insert(insertion_pos, score)
            self.items.insert(insertion_pos, item)
            if len(self.items) > self.num_to_keep:
                self.items = self.items[:-1]
                self.scores = self.scores[:-1]

    def get_best(self):
        return self.scores, self.items



