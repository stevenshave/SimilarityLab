"""
Generate MRLogP descriptors from input SMILES

The MRLogP package requires molecules represented by 3 sets of molecular
descriptors.
This requires OpenBabel The easiest way to do this is with:
'conda install -c conda-forge openbabel'
By Steven Shave - stevenshave@gmail.com
"""
from openbabel import openbabel
import argparse
from rdkit import Chem
from pathlib import Path
from rdkit import RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRCAT
from utils.rdkonf6 import smiles_to_3dmol
import sys, os, pickle
import numpy as np
import pandas as pd
import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
with tf.device("cpu:0"):
   from tensorflow.compat.v1.keras.backend import set_session
   from tensorflow.compat.v1.keras import backend as K
   from tensorflow.compat.v1.keras.models import Sequential
   from tensorflow.compat.v1.keras.layers import Dense
   from tensorflow.compat.v1.keras.layers import Dropout
   from tensorflow.compat.v1.keras.optimizers import Adam
   from tensorflow.compat.v1.keras.callbacks import Callback
   from tensorflow.compat.v1.keras.callbacks import ModelCheckpoint
   from tensorflow.compat.v1.keras.layers import PReLU
   from sklearn.preprocessing import StandardScaler
   from tensorflow.compat.v1.keras.models import load_model

class MRLogPPredictor():
    
    def __init__(self, scaler_pickle='/home/ubuntu/similarity_lab/utils/scaler.pkl'):
        config = tf.ConfigProto()
        set_session(tf.Session(config=config))
        
        self.scaler = pickle.load(open(scaler_pickle,'rb'))

    def predict_from_smiles(self, smiles, model_hdf5="/home/ubuntu/similarity_lab/utils/model-tl-18.hdf5"):
        if " " not in smiles:
            smiles=smiles+" querymol"
        smiles, title=smiles.split()
        descriptors=self.get_mrlogP_descriptors(smiles, title)
        def root_mean_squared_error(y_true, y_pred):
            return K.sqrt(K.mean(K.square(y_pred - y_true)))
        self.model = load_model(model_hdf5, custom_objects={'root_mean_squared_error': root_mean_squared_error})
        return self.predict_from_descriptors(descriptors)
    
    def predict_from_descriptors(self,mrlogp_descriptors):
      mrlogp_descriptors[256:]=self.scaler.transform(mrlogp_descriptors[256:].reshape(1,-1))[0]
      return self.model.predict(mrlogp_descriptors.reshape(1,-1))[0][0]

    def get_mrlogP_descriptors(self, smiles:str, moltitle:str):
        # Descriptors are morgan, fp4, then USRCAT
        mrlogP_descriptor_length=128+128+60
        
        obConversion=openbabel.OBConversion()
        obConversion.SetInAndOutFormats("smi", "mdl")
        ob_mol=openbabel.OBMol()

        # Create RDKit and OpenBabel molecules
        rdkit_mol=Chem.AddHs(smiles_to_3dmol(smiles,"querymol"))
        obConversion.ReadString(ob_mol, smiles)
        
        # Generate Morgan/ECFP4
        morgan_fingerprint = AllChem.GetMorganFingerprintAsBitVect(Chem.RemoveHs(rdkit_mol),2,128).ToBitString()
        # Generate USRCAT
        usrcat_descriptors= GetUSRCAT(rdkit_mol)
        # Generate FP4
        fp4fp = openbabel.vectorUnsignedInt()
        fingerprinter = openbabel.OBFingerprint.FindFingerprint("FP4")
        
        fingerprinter.GetFingerprint(ob_mol, fp4fp)

        openbabel.OBFingerprint.Fold(fingerprinter,fp4fp, 128)

        logP_descriptors=np.full((mrlogP_descriptor_length), np.nan)
        
        for i,v in enumerate(morgan_fingerprint):
            logP_descriptors[i]=float(v)
        
        fp4_p1=[float(x) for x in list(format(fp4fp[0],'032b'))]
        fp4_p2=[float(x) for x in list(format(fp4fp[1],'032b'))]
        fp4_p3=[float(x) for x in list(format(fp4fp[2],'032b'))]
        fp4_p4=[float(x) for x in list(format(fp4fp[3],'032b'))]
        logP_descriptors[128:256]=fp4_p1+fp4_p2+fp4_p3+fp4_p4

        for i,v in enumerate(usrcat_descriptors):
            logP_descriptors[256+i]=float(v)
        return logP_descriptors

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("smiles", help="Smiles file")
    args = parser.parse_args()
    mrlogP=MRLogPPredictor()
    print("SMILES", args.smiles)
    print(mrlogP.predict_from_smiles(args.smiles))