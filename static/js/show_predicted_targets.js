function createPredictedTargetsTable() {

    tableHeadings = ["Target", "Hit count", "ChEMBL IDs of actives"]
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
            let text = document.createTextNode(d);
                cell.appendChild(text);
            
        });
    });
}
