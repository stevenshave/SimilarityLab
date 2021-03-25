import os
from pathlib import Path

class Config(object):
    SECRET_KEY = os.environ.get('FLASK_SECRET_KEY') or '3Sy63TSRyVBaq6G37Ch7lXYMLQPuYgkUfRAN1Av4ltc'
    CELERY_BROKER_URL = 'redis://localhost:6379/0'
    #CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
    DATASETS_DIRECTORY=Path("/home/ubuntu/data/")
    DATASETS=[
        [99, "10ktestset", "Small set of 10k druglike molecules (QED scores > 0.9)"],
        [1, "eMolecules202103", "eMolecules-2021-03 (29,368,630 unique molecules)"],

    ]
    QUERY_SIMILARS_DIRECTORY=DATASETS_DIRECTORY/"queries/similars/"
    QUERY_TARGETS_DIRECTORY=DATASETS_DIRECTORY/"queries/targets/"

    CHEMBL_VERSION_NUMBER=28
    CHEMBL_USRCATSL_ROOT=DATASETS_DIRECTORY/"chembl28_actives.sdf.usrcatsl"
    CHEMBL_USRCATSL_BIN=Path(str(CHEMBL_USRCATSL_ROOT)+".bin")
    assert(CHEMBL_USRCATSL_BIN.exists()),"Missing bin file"
    CHEMBL_USRCATSL_SMI=Path(str(CHEMBL_USRCATSL_ROOT)+".smi")
    assert (CHEMBL_USRCATSL_SMI), "Missing smiles file"
    NUM_TO_KEEP_SIMILARS=200
    NUM_TO_KEEP_TARGETS=100

    CCHEMBLID_TO_TCHEMBLIDS_PATH=DATASETS_DIRECTORY/"cchemblid_to_tchemblids.json"
    assert CCHEMBLID_TO_TCHEMBLIDS_PATH.exists(), "Missing cchemblid_to_tchemblids.json"
    TCHEMBLIDS_TO_PREFNAMES_PATH=DATASETS_DIRECTORY/"tchemblids_to_prefnames.json"
    assert TCHEMBLIDS_TO_PREFNAMES_PATH.exists(), "Missing tchemblids_to_prefnames.json"
    
    CCHEMBLID_TO_TCHEMBLIDS=None
    TCHEMBLIDS_TO_PREFNAMES=None
