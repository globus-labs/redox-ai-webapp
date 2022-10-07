'use strict';

function compute() {
    // Get the name of the molecule
    var mol = document.getElementById("smiles").value;

    axios.get("/api/render", {params: {smiles: mol}})
        .then(function (res) {
            document.getElementById("mol-img").innerHTML = res.data
        })
}