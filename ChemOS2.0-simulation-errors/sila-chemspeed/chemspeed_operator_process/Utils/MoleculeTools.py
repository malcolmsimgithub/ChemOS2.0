from openbabel import openbabel
obabel_converter = openbabel.OBConversion()
obabel_converter.SetInAndOutFormats("smi", "mdl")


def formula_from_smiles(smiles):
    mol = openbabel.OBMol()
    obabel_converter.ReadString(mol, smiles)
    return mol.GetFormula()


def exact_mass_from_smiles(smiles):
    mol = openbabel.OBMol()
    obabel_converter.ReadString(mol, smiles)
    return mol.GetExactMass()


def mw_from_smiles(smiles):
    mol = openbabel.OBMol()
    obabel_converter.ReadString(mol, smiles)
    return mol.GetMolWt()

# TODO: Move to deprecated tools after restructuring the architecture -> not needed any more
