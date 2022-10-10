'use strict';

// We use a websocket to send results from the web server as they complete
var socket = new WebSocket("ws://127.0.0.1:8000/ws")

// Get the current molecule
function activeMolecule() {
    return document.getElementById("smiles").value;
}

// Function that begins a computation
function compute() {
    // Get the name of the molecule
    var mol = activeMolecule();

    // Render it as an SVG
    axios.get("/api/render", {params: {smiles: mol}})
        .then(function (res) {
            document.getElementById("mol-img").innerHTML = res.data
        })

    // Request its properties
    socket.send(JSON.stringify({"smiles": mol}))
}


// Functions
socket.onmessage = function(event) {
    /** Set the appropriate molecule property **/
    var data = JSON.parse(event.data);
    var smiles = data.smiles;

    // Make sure the molecule hasn't changed
    var set_smiles = activeMolecule();
    if (set_smiles != smiles) {
        console.log('Molecule has updated')
        return
    }

    // Update the appropriate row
    var model_id = data.model_name;
    var value = data.value;
    var result_row = document.getElementById(model_id);
    if (! result_row) {
        return console.log(`[error] result row not found ${model_id}`)
    }
    result_row.cells[1].innerHTML = value
}

