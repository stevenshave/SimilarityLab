from utils.cheminformatics_tools_std_include import *
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

usrcat_binary_file=open("/home/ubuntu/data/chembl28_actives.sdf.usrcatsl.bin", "rb")
struct_size_in_bytes=struct.calcsize(usrcat_binary_struct_format_string)
usrcat_binary_file.seek(0, io.SEEK_END)
binary_file_size = usrcat_binary_file.tell()
usrcat_binary_file.seek(0)
num_mols_in_binary_file=binary_file_size//struct_size_in_bytes

query_descriptors=[0.0]*60
print(query_descriptors)
keepn = KeepNHighest(200)
for molcounter in range(num_mols_in_binary_file):
    pos_and_desc_bytes=usrcat_binary_file.read(struct_size_in_bytes)
    descriptors=struct.unpack_from(usrcat_binary_struct_format_string, pos_and_desc_bytes,0)
    smiles_line_number=descriptors[0]
    usrcat_score = GetUSRScore(query_descriptors, descriptors[1:])
    keepn.insert(smiles_line_number, usrcat_score)
best_scores, best_smiles_line_numbers = keepn.get_best()
linecounter=1
best_smiles_line_numbers_set=set(best_smiles_line_numbers)
smiles_dict={}
with open("/home/ubuntu/data/chembl28_actives.sdf.usrcatsl.smi") as smiles_file:
    line=smiles_file.readline()
    while line:
        if linecounter in best_smiles_line_numbers_set:
            smiles_dict[linecounter]=line.rstrip()
        line=smiles_file.readline()
        linecounter+=1

print("Done")