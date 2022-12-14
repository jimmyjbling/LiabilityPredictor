from rdkit.Chem import MolFromSmiles, Draw
from LiabilityPredictor.assay_liability_calculator import main

COLORS = {
    0: "green",
    1: "#b09b2a",
    2: "red",
    "Possible Interference": "#b09b2a",
    "No Interference": "green",
    "Putative aggregator": "#b09b2a",
    "Non-aggregator": "green"
}


def get_molecule_data_from_smiles(smiles_str, options):
    molecule = MolFromSmiles(smiles_str)

    if molecule is None:
        return None

    data = main(smiles_str, **options)

    # add color coding of text
    data = [_ + [COLORS[_[1]]] for _ in data]

    # skip if no models selected
    if len(data) == 0:
        return None

    return {
        'svg': Draw.MolsToGridImage([molecule], useSVG=True, molsPerRow=1),
        'SMILES': smiles_str,
        'pred_data': data,
    }
