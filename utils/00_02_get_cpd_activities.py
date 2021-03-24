"""Make cpd_chemblid to targets hit chemblid lookup table

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
from progressbar import progressbar
from pathlib import Path
from easy_chembl_queries import *
import json
output_file=open("actives.smi", "wt")

cpdchemblids=list(set([cid.strip() for line in progressbar(open("/home/ubuntu/data/chembl28_actives.sdf.usrcatsl.smi").readlines()) for cid in line.split(" ")[1].split(";")]))

chembl = chemblquerier.ChEMBLQuerier(Path("/home/ubuntu/data/chembl_28.db"))

cpdchemblid_to_molregno={}
for cpdchemblid in progressbar(cpdchemblids):
    molregno=chembl.query("select molregno from molecule_dictionary where chembl_id==\""+cpdchemblid+"\"")[0][0]
    cpdchemblid_to_molregno[cpdchemblid]=molregno

molregno_to_assayids={}
for molregno in progressbar(cpdchemblid_to_molregno.values()):
    assay_ids=[aid[0] for aid in chembl.query("select distinct(assay_id) from activities where molregno=="+str(molregno)+" and standard_value<=10000 and standard_type in ('IC50', 'Kd', 'MIC', 'EC50', 'ED50', 'KD')") if len(aid)>0]
    molregno_to_assayids[molregno]=assay_ids

assay_ids_to_tids={}
for assay_id in progressbar(set([aid for assay_id_lists in molregno_to_assayids.values() for aid in assay_id_lists])):
    tid=chembl.query("select tid from assays where assay_id="+str(assay_id))[0][0]
    assay_ids_to_tids[assay_id]=tid

tids_to_tchemblids={}
tchemblids_to_prefnames={}

for tid in set([tid for tid in assay_ids_to_tids.values()]):
    cid_and_prefname_list=(chembl.query("select chembl_id, pref_name from target_dictionary where tid="+str(tid)+" and target_type in ('SINGLE PROTEIN','PROTEIN','PROTEIN-PROTEIN INTERACTION','PROTEIN FAMILY')"))
    if len(cid_and_prefname_list)==0:continue
    tids_to_tchemblids[tid]=cid_and_prefname_list[0][0]
    tchemblids_to_prefnames[cid_and_prefname_list[0][0]]=cid_and_prefname_list[0][1]

cchemblid_to_tchemblids={}
for cchemblid in cpdchemblid_to_molregno.keys():
    cchemblid_to_tchemblids[cchemblid]=set()
    molregno=cpdchemblid_to_molregno[cchemblid]
    for assay_id in molregno_to_assayids[molregno]:
        tid=assay_ids_to_tids[assay_id]
        if tid in tids_to_tchemblids.keys():
            cchemblid_to_tchemblids[cchemblid].add(tids_to_tchemblids[tid])
    cchemblid_to_tchemblids[cchemblid]=list(cchemblid_to_tchemblids[cchemblid])

json.dump(cchemblid_to_tchemblids,open("/home/ubuntu/data/cchemblid_to_tchemblids.json", "w"))
json.dump(tchemblids_to_prefnames,open("/home/ubuntu/data/tchemblids_to_prefnames.json", "w"))


exit()
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
