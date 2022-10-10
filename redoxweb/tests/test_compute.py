from redoxweb.compute import LogPModel, MolecularWeightModel


def test_logp():
    model = LogPModel()
    assert model.run("C").value == 0.6361


def test_molwt():
    model = MolecularWeightModel()
    assert model.run("C").value == 16.043
