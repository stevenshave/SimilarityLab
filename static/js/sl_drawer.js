var composer = new Kekule.Editor.Composer(document.getElementById('molecule-drawer'))
composer.setPredefinedSetting('molOnly');

function get_smiles() {
    var smiles_list=[];
    var molecules=composer.exportObjs(Kekule.Molecule);
    var i;
    for (i = 0; i < molecules.length; i++) {
        smiles_list.push(Kekule.IO.saveFormatData(molecules[i],"smi"));
    }
    
    smiles=smiles_list.join(".");
    document.getElementById("molsmiles").value=smiles;
}

function get3Dsimilars(){
    databaseChoice=document.getElementById("databaseChoice").value;
    window.location.href = "show_similars.html";
}
function predictTargets(){
    window.location.href = "show_predicted_targets.html";
}
function exploreCompound(){
    window.location.href = "show_compound.html";
}

composer.on('editObjsUpdated', function(){
    get_smiles()
  });