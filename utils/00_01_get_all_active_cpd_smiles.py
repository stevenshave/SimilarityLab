"""Get SMILES of actives

Get all actives in ChEMBL - defin4ed by an activity to sigular protein target of at least 10 µM IC50/KD

'IC50', 'Kd', 'MIC', 'EC50', 'ED50', 'KD'
As of ChEMBL28- the following activity counts are present for activities with standard_values<=10000:
(1483127, 'IC50'),
(1472916, 'Potency'),
(1024210, 'Inhibition'),
(586019, 'GI50'),
(571403, 'Activity'),
(557449, 'MIC'),
(540653, 'Ki'),
(339133, 'INHIBITION'), 
(286524, 'EC50'),
(104002, 'Ratio IC50'),
(100661, 'Ratio'),
(94587, 'ED50'),
(91258, 'GI'),
(88244, 'IZ'), 
(86266, 'Kd'),
(73944, 'Residual Activity'),
(62207, 'T1/2'),
(61452, 'AC50'), 
(51736, 'CL'), 
(44460, '% Control'),
(43001, 'FC'),
(35371, 'Ratio Ki')
(30813, 'Drug uptake')
(29424, 'MIC90'), 
      


"""
import progressbar
from pathlib import Path
from easy_chembl_queries import *
import json
output_file=open("actives.smi", "wt")

chembl = chemblquerier.ChEMBLQuerier(Path("d:\chembl_28.db"))
single_protein_tids = [f"'{x[0]}'" for x in chembl.query(
    "select tid from target_dictionary where target_type in ('SINGLE PROTEIN','PROTEIN','PROTEIN-PROTEIN INTERACTION','PROTEIN FAMILY')"
)]
assay_ids=[f"'{x[0]}'" for x in chembl.query("select assay_id from assays where tid in ("+",".join(single_protein_tids)+")")]
molregs = [x[0] for x in chembl.query("select distinct(molregno) from activities where standard_value<=10000 and standard_type in ('IC50', 'Kd', 'MIC', 'EC50', 'ED50', 'KD') and assay_id in ("+",".join(assay_ids)+")")]
for molregno in progressbar.progressbar(molregs):
    chemblid=chembl.query("select chembl_id from molecule_dictionary where molregno='"+str(molregno)+"'")
    if len(chemblid)==0:continue
    chemblid=chemblid[0]
    if len(chemblid)==0:continue
    smiles=chembl.query("select canonical_smiles from compound_structures where molregno='"+str(molregno)+"'")
    if len(smiles)==0:continue
    smiles=smiles[0]
    if len(smiles)==0:continue

    output_file.write(f"{smiles[0]} {chemblid[0]}\n")
output_file.close()
# print("Getting active human protein chemblids")
# human_protein_chemblids=[p[0] for p in chembl.get_human_proteins()]
# print("Got", len(human_protein_chemblids), "human protein ids")

# molregnumbers=chembl.get_cpds_molreg_active_against_chemblid_target_list([str(cid) for cid in human_protein_chemblids])
# smiles_molreg_list=chembl.get_smiles_for_molregs(molregnumbers)

# json.dump(smiles_molreg_list, open("dat_active_smiles_with_molreg.json", "w"))
