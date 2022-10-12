from redoxweb.compute import LogPModel, MolecularWeightModel, SolvationEnergyModel


def test_logp():
    model = LogPModel()
    assert model.run("C").value == 0.6361


def test_molwt():
    model = MolecularWeightModel()
    assert model.run("C").value == 16.043


def test_solvation():
    model = SolvationEnergyModel()
    model.json()

    assert model.run("C").value_str == "0.51"
