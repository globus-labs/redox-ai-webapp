"""Confiuration for the web service"""
from redoxweb.compute import LogPModel, MolecularWeightModel

# TODO (wardlt): Make this read from a YAML file

models = [MolecularWeightModel(), LogPModel()]
