function createResultsTable() {

    tableHeadings = ["2D/SMILES", "USRCAT Score", "Morgan score", "eMolecules ID", "MW"]
    myTable = document.getElementById("resultsTable")

    //Generate the head
    let thead = myTable.createTHead();
    let row = thead.insertRow();
    tableHeadings.forEach(function (heading) {
        let th = document.createElement("th");
        let text = document.createTextNode(heading);
        th.appendChild(text);
        row.appendChild(th);
    });

    //Generate and fill the table data

    moleculeData.forEach(function (moldata, moliterator) {
        let row = myTable.insertRow();
        moldata.forEach(function (d, di) {
            let cell = row.insertCell();
            if (di == 0) {
                cell.innerHTML = "<canvas id='canvas" + moliterator + "'></canvas><br>" + moldata[0];
            } else {
                let text = document.createTextNode(d);
                cell.appendChild(text);
            }
        });
    });
}
function drawMolecules() {
    moleculeData.forEach(function (moldata, moliterator) {
        SmilesDrawer.parse(moldata[0], function (tree) {
            smilesDrawer.draw(tree, 'canvas' + moliterator, 'light', false);
        }, function (err) {
            console.log(err);
        });
        canv = document.getElementById('canvas' + moliterator);

    });


}