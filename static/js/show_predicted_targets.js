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
        let cell0 = row.insertCell();
        let cell1 = row.insertCell();
        let cell2 = row.insertCell();
        cell0.innerHTML="<a href='https://www.ebi.ac.uk/chembl/g/#search_results/targets/query="+moldata[0]+"'>"+moldata[0]+"</a>";
        cell1.innerHTML=moldata[1];
        cell2.innerHTML ="";
        chemblids=moldata[2];
        
        const myArr = chemblids.split(", ");
        for (var i = 0; i < myArr.length; i++)
            cell2.innerHTML = cell2.innerHTML+"<a href='https://www.ebi.ac.uk/chembl/compound_report_card/"+myArr[i]+"'>"+myArr[i]+"</a>, ";
        
    });
}
