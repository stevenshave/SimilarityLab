import os, io, datetime, struct
from pathlib import Path
from rdkit.Chem.inchi import MolToInchiKey
from utils.cheminformatics_tools_std_include import standardise_mol_remove_salts, standardise_smiles_remove_salts, usrcat_binary_struct_format_string, KeepNHighest
from rdkit.Chem.AllChem import GetMorganFingerprint
from rdkit.DataStructs import DiceSimilarity
from flask.globals import request
from config import Config
from celery import Celery
from flask import Flask, render_template, url_for, send_from_directory, redirect, send_file
from forms import *
from rdkit import Chem

from utils.rdkonf6 import smiles_to_3dmol
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT
from rdkit.Chem.Descriptors import MolWt

app=Flask(__name__)
app.config.from_object(Config)
print(app.config['SECRET_KEY'])


# Celery queue client setup
celery_client = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery_client.conf.update(app.config)

@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.route('/')
@app.route('/index.html')
def index():
    return render_template("index.html")


@app.route("/download_similars.csv", methods=["GET"])
def download_similars():
    mol_inchi_and_dbid=request.args.get('molid')
    
    print("Here")
    # Check people are not trying to access other files
    if len(mol_inchi_and_dbid.split("_")[0])!=27 or "/" in mol_inchi_and_dbid or "." in mol_inchi_and_dbid or "~" in mol_inchi_and_dbid:
        return redirect("/")
    
    return send_file(os.path.join(app.config['QUERY_SIMILARS_DIRECTORY'], mol_inchi_and_dbid+".csv"), as_attachment=True,attachment_filename="similars.csv")
    #return send_from_directory(app.config['QUERY_SIMILARS_DIRECTORY'],
    #                           mol_inchi_and_dbid+".csv", mimetype='text/csv', attachment_filename="similars.csv")


@app.route("/show_similars.html", methods=["GET"])
def show_similars():
    
    #Parse the get request for the mol argument
    mol_inchi_and_dbid=request.args.get('mol')
    
    # No molecule inchi and DB ID supplied
    if mol_inchi_and_dbid is None:
        return redirect("/")

    # Inchi and DB supplied- check if 1) Submitted? 2) Finished?
    if (Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(mol_inchi_and_dbid+".info")).exists():
        if not (Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(mol_inchi_and_dbid+".csv")).exists():
            return render_template("message.html", heading="Not yet complete", message="Job is submitted, please check back later.")
        # Must be complete, CSV exists, extract lines and pass to rendering of show_similars
        # Has  to look like this:
        #    ["CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", 0.9, 0.6,"12345", "206.285"],
        #    ["O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl", 0.82, 0.4, "54321", "296.15"],
        #    ["O=P(O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3O", 0.3,0.9, "54321", "347.22"],
        #
        # Where the above is SMILES, USRCATScore, MorganScore, ID, Molecular weight

        moldat=""
        for line in open(Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(mol_inchi_and_dbid+".csv")).readlines()[1:]:
            smiles, usrcatscore, morganscore, id, mw=line.split(",")
            moldat+=f"[\"{smiles}\",{usrcatscore},{morganscore},\"{id}\",{mw}],\n"
        return render_template("show_similars.html", moldata=moldat, mol_inchi_and_dbid=mol_inchi_and_dbid )
    else:
        return render_template("message.html", heading="Not found", message="It seems that job was never submitted.")


@app.route('/find_similars.html',  methods=['GET', 'POST'])
def find_similars():
    form=FindSimilarsForm()
    form.select.choices=[(ds[0], ds[2]) for ds in app.config["DATASETS"]]

    if form.validate_on_submit():
        # Correct IP courtesy of https://stackoverflow.com/questions/3759981/get-ip-address-of-visitors-using-flask-for-python
        # Get requesting IP:
        client_ip=None
        if request.environ.get('HTTP_X_FORWARDED_FOR') is None:client_ip=request.environ['REMOTE_ADDR']
        else:client_ip=request.environ['HTTP_X_FORWARDED_FOR'] # if behind a proxy
        smiles_std=standardise_smiles_remove_salts(form.smiles.data)
        mol=smiles_to_3dmol(smiles_std)
        if mol.GetNumHeavyAtoms()<3:
            return render_template("message.html", heading="Molecule too small", message="The USRCAT molecular similarity technique requires molecules to be composed of at least 3 heavy atoms.")
            
        inchi_key=Chem.inchi.MolToInchiKey(mol)
        query_mol_identifier=inchi_key+"_"+form.select.data
        if mol is None:
            return render_template("message.html", heading="3D generation info", message="3D generation failed. This could be beacuse it is too big, too flexible, the submitted SMILES is invalid, or contains metals that SimilarityLab (using the 3D generation method detailed in the about section) is unable to find parameters for. Allowed atom types are: C, N, O, S, F, Cl, Br, I, B, P, Si, and H.")
        usrcat_descriptors=GetUSRCAT(mol)
    
        database_binary_path=app.config['DATASETS_DIRECTORY']/Path(app.config['DATASETS'][int(form.select.data)][1]+".sdf.usrcatsl.bin")
        database_smiles_path=app.config['DATASETS_DIRECTORY']/Path(app.config['DATASETS'][int(form.select.data)][1]+".sdf.usrcatsl.smi")
        
        # Check if a cached version exists before firing off a new request. If it does, then just show it.
        if (Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(query_mol_identifier+".info")).exists():
            with open(Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(query_mol_identifier+".info"),"a") as infofile:
                infofile.write(f"{query_mol_identifier},{smiles_std},{client_ip},{datetime.datetime.now()}\n")
            return redirect(url_for("show_similars")+"?mol="+query_mol_identifier)
        with open(Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(query_mol_identifier+".info"), "w") as infofile:
            infofile.write(f"{query_mol_identifier},{smiles_std},{client_ip},{datetime.datetime.now()}\n")
        # Not found, so we make a request
        print("Going to call GET SIMILAR MOLECULES")
        get_similar_molecules.apply_async(args=[usrcat_descriptors,  smiles_std, str(database_binary_path), str(database_smiles_path), int(form.select.data)], countdown=1)
        return render_template("message.html", heading="Finding similars", message="The database is currently being queried for similars to your uploaded molecule ("+smiles_std+").<br>Please check the link bellow periodically to view your results. Searches against the entire ~26M eMolecules can take up to 2 hrs and even longer when the server is under heavy load. <br><a href='"+url_for("show_similars")+"?mol="+inchi_key+"_"+form.select.data+"'>Click here to check status</a>")

    for fieldName, errorMessages in form.errors.items():
        for err in errorMessages:
            print(err, fieldName)
    return render_template("find_similars.html", form=form)

@app.route('/predict_targets.html',  methods=['GET', 'POST'])
def predict_targets():
    form=PredictTargetsForm()

    if form.validate_on_submit():
        # Correct IP courtesy of https://stackoverflow.com/questions/3759981/get-ip-address-of-visitors-using-flask-for-python
        # Get requesting IP:
        client_ip=None
        if request.environ.get('HTTP_X_FORWARDED_FOR') is None:client_ip=request.environ['REMOTE_ADDR']
        else:client_ip=request.environ['HTTP_X_FORWARDED_FOR'] # if behind a proxy
        smiles_std=standardise_smiles_remove_salts(form.smiles.data)
        mol=smiles_to_3dmol(smiles_std)
        if mol is None:
            return render_template("message.html", heading="3D generation info", message="3D generation failed. This could be beacuse it is too big, too flexible, the submitted SMILES is invalid, or contains metals that SimilarityLab (using the 3D generation method detailed in the about section) is unable to find parameters for. Allowed atom types are: C, N, O, S, F, Cl, Br, I, B, P, Si, and H.")
        if mol.GetNumHeavyAtoms()<3:
            return render_template("message.html", heading="Molecule too small", message="Molecule is too small, please query at least 3 heavy atoms.")
            
        inchi_key=Chem.inchi.MolToInchiKey(mol)
        usrcat_descriptors=GetUSRCAT(mol)
    
        # Check if a cached version exists before firing off a new request. If it does, then just show it.
        if (Path(app.config['QUERY_TARGETS_DIRECTORY'])/(inchi_key+".info")).exists():
            with open(Path(app.config['QUERY_TARGETS_DIRECTORY'])/(inchi_key+".info"),"a") as infofile:
                infofile.write(f"{inchi_key},{smiles_std},{client_ip},{datetime.datetime.now()}\n")
            return redirect(url_for("show_predicted_targets.html")+"?mol="+inchi_key)
        with open(Path(app.config['QUERY_TARGETS_DIRECTORY'])/(inchi_key+".info"), "w") as infofile:
            infofile.write(f"{inchi_key},{smiles_std},{client_ip},{datetime.datetime.now()}\n")
        # Not found, so we make a request
        get_predicted_targets.apply_async(args=[usrcat_descriptors,  smiles_std, str(app.config['CHEMBL_USRCATSL_BIN']), str(app.config['CHEMBL_USRCATSL_SMI'])], countdown=1)
        return render_template("message.html", heading="Finding similars", message="ChEMBL is being queried for similars to your uploaded molecule ("+smiles_std+").<br>Please check the link bellow periodically to view your results. Searches against the entire ~26M eMolecules can take up to 2 hrs and even longer when the server is under heavy load. <br><a href='"+url_for("show_similars")+"?mol="+inchi_key+"_"+form.select.data+"'>Click here to check status</a>")

    for fieldName, errorMessages in form.errors.items():
        for err in errorMessages:
            print(err, fieldName)
    return render_template("predict_targets.html", form=form)




@app.route('/explore_compound.html')
def explore_compound():
    return render_template("explore_compound.html")




from subprocess import Popen, PIPE



import time
@celery_client.task
def get_similar_molecules(query_descriptors:list, query_smiles:str, database_binary_path:str, database_smiles_path:str, database_id:int):
    """Celery task that reads database binary files comparing query descriptors

    Args:
        query_descriptors (list): USRCAT descriptors of query.
        query_smiles (str): SMILES representation of the query molecule.
        database_binary_path (str): Binary file path - must be string, not path, as Path is non-serialisable.
        database_smiles_path (str): Smiles file path - same as above, must be a string.
        database_id (int): Int representing database ID, used to cache results against specific databases
    """    
    print("Worker running for "+ query_smiles)
    mol=Chem.MolFromSmiles(query_smiles)
    query_mol_identifier=Chem.inchi.MolToInchiKey(mol)+"_"+str(database_id)


    # CPP program bellow called for speed of processing
    ###################################################################################
    # Build the command line -    //Arguments must be 
    # 0: Executable
    # 1: Binary file location without last .bin extension, so that .bin and .smi file locations can be derived
    # 2: Number of best to keep
    # 3-63: USRCAT descriptors of query
    command_line=["/home/ubuntu/similarity_lab/utils/usrcat_binary_reader_similarity_lab", database_binary_path.replace(".bin",""), str(app.config['NUM_TO_KEEP'])]
    for i in range(60):
        command_line.append(str(query_descriptors[i]))
    process = Popen(command_line, stdout=PIPE)
    output, err = process.communicate()
    exit_code = process.wait()
    lines=output.decode("utf-8").splitlines()

    with open(Path(app.config['QUERY_SIMILARS_DIRECTORY'])/(query_mol_identifier+".csv"),"w") as outfile:
        for line in lines:
            stripped_line=line.strip()
            candidate_smiles, title_comma_score=stripped_line.split(" ")
            candidate_title, candidate_score = title_comma_score.split(",")
            candidate_score=float(candidate_score)
            candidate_mol=Chem.MolFromSmiles(candidate_smiles)
            candidate_morgan_score=DiceSimilarity(GetMorganFingerprint(candidate_mol,2),GetMorganFingerprint(mol,2))
            mw=MolWt(candidate_mol)
            outfile.write(f'{candidate_smiles},{round(candidate_score,3)},{round(candidate_morgan_score,3)},{candidate_title.replace("_1","")},{round(mw,3)}\n')
    print("Worker done")

@celery_client.task
def get_predicted_targets(query_descriptors:list, query_smiles:str, database_binary_path:str, database_smiles_path:str):
    """Celery task that reads database binary files comparing query descriptors

    Args:
        query_descriptors (list): USRCAT descriptors of query.
        query_smiles (str): SMILES representation of the query molecule.
        database_binary_path (str): Binary file path - must be string, not path, as Path is non-serialisable.
        database_smiles_path (str): Smiles file path - same as above, must be a string.
        database_id (int): Int representing database ID, used to cache results against specific databases
    """    
    
    print("Worker running for "+ query_smiles)
    mol=Chem.MolFromSmiles(query_smiles)
    query_mol_identifier=Chem.inchi.MolToInchiKey(mol)+"_"+str(app.config['CHEMBL_VERSION_NUMBER'])


    # CPP program bellow called for speed of processing
    ###################################################################################
    # Build the command line -    //Arguments must be 
    # 0: Executable
    # 1: Binary file location without last .bin extension, so that .bin and .smi file locations can be derived
    # 2: Number of best to keep
    # 3-63: USRCAT descriptors of query
    command_line=["/home/ubuntu/similarity_lab/utils/usrcat_binary_reader_similarity_lab", database_binary_path.replace(".bin",""), str(app.config['NUM_TO_KEEP'])]
    for i in range(60):
        command_line.append(str(query_descriptors[i]))
    process = Popen(command_line, stdout=PIPE)
    output, err = process.communicate()
    exit_code = process.wait()
    lines=output.decode("utf-8").splitlines()

    with open(Path(app.config['QUERY_TARGETS_DIRECTORY'])/(query_mol_identifier+".csv"),"w") as outfile:
        for line in lines:
            print(line)
            stripped_line=line.strip()
            candidate_smiles, title_comma_score=stripped_line.split(" ")
            candidate_title, candidate_score = title_comma_score.split(",")
            candidate_score=float(candidate_score)
            candidate_mol=Chem.MolFromSmiles(candidate_smiles)
            candidate_morgan_score=DiceSimilarity(GetMorganFingerprint(candidate_mol,2),GetMorganFingerprint(mol,2))
            mw=MolWt(candidate_mol)
            outfile.write(f'{candidate_smiles},{round(candidate_score,3)},{round(candidate_morgan_score,3)},{candidate_title.replace("_1","")},{round(mw,3)}\n')
    print("Worker done")




if __name__ == "__main__":
    app.run(debug=True)
