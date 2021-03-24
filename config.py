import os
from pathlib import Path

class Config(object):
    SECRET_KEY = os.environ.get('FLASK_SECRET_KEY') or '3Sy63TSRyVBaq6G37Ch7lXYMLQPuYgkUfRAN1Av4ltc'
    CELERY_BROKER_URL = 'redis://localhost:6379/0'
    #CELERY_RESULT_BACKEND = 'redis://localhost:6379/0'
    DATASETS_DIRECTORY=Path("/home/ubuntu/data/")
    DATASETS=[
        [0, "10ktestset", "Small set of 10k druglike molecules (QED scores > 0.9)"],
    ]
    QUERY_SIMILARS_DIRECTORY=Path("/home/ubuntu/data/queries/similars/")
    QUERY_TARGETS_DIRECTORY=Path("/home/ubuntu/data/queries/targets/")

    CHEMBL_VERSION_NUMBER=28
    CHEMBL_USRCATSL_ROOT=Path("/home/ubuntu/data/chembl28_actives.sdf.usrcatsl")
    CHEMBL_USRCATSL_BIN=Path(str(CHEMBL_USRCATSL_ROOT)+".bin")
    CHEMBL_USRCATSL_SMI=Path(str(CHEMBL_USRCATSL_ROOT)+".smi")

    NUM_TO_KEEP=200
