from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import SimilarityMaps
from scipy.spatial.distance import cdist
import numpy as np

import glob
import os
import pickle
import io
import matplotlib.pyplot as plt

# import warnings
# def warn(*args, **kwargs):
#     pass
# warnings.warn = warn


MODEL_DICT = {
    'Firefly_Luciferase_counter_assay_training_set_curated_model.pkl': 'Firefly Luciferase interference',
    'Nano_Luciferase_counter_assay_training_set_curated_model.pkl': 'Nano Luciferase interference',
    'Redox_training_set_curated_model.pkl': 'Redox interference',
    'Thiol_training_set_curated_model.pkl': 'Thiol interference',
    'agg_betalac_model.pkl': 'AmpC β-lactamase aggregation',
    'agg_cruzain_model.pkl': 'Cysteine protease cruzain aggregation'
}

MODEL_DICT_INVERT = {val: key for key, val in MODEL_DICT.items()}

OUTCOME_DICT = {
    'Firefly Luciferase interference': {
        1: "Possible Interference",
        0: "No Interference"
    },
    'Nano Luciferase interference': {
        1: "Possible Interference",
        0: "No Interference"
    },
    'Redox interference': {
        1: "Possible Interference",
        0: "No Interference"
    },
    'Thiol interference': {
        1: "Possible Interference",
        0: "No Interference"
    },
    'AmpC β-lactamase aggregation': {
        1: "Putative aggregator",
        0: "Non-aggregator"
    },
    'Cysteine protease cruzain aggregation': {
        1: "Putative aggregator",
        0: "Non-aggregator"
    }
}


def run_prediction(model, model_data, smiles, calculate_ad=True):
    fp = np.zeros((2048, 1))
    _fp = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smiles), radius=3, nBits=2048)
    DataStructs.ConvertToNumpyArray(_fp, fp)

    pred_proba = model.predict_proba(fp.reshape(1, -1))[:, 1]
    pred = 1 if pred_proba > 0.5 else 0

    if pred == 0:
        pred_proba = 1-pred_proba

    if calculate_ad:
        ad = model_data["D_cutoff"] > np.min(cdist(model_data['Descriptors'].to_numpy(), fp.reshape(1, -1)))
        return pred, pred_proba, ad
    return pred, pred_proba, None


def get_prob_map(model, smiles):
    def get_fp(mol, idx):
        fps = np.zeros((2048, 1))
        _fps = SimilarityMaps.GetMorganFingerprint(mol, idx, radius=3, nBits=2048)
        DataStructs.ConvertToNumpyArray(_fps, fps)
        return fps

    def get_proba(fps):
        return float(model.predict_proba(fps.reshape(1, -1))[:, 1])

    mol = Chem.MolFromSmiles(smiles)
    fig, _ = SimilarityMaps.GetSimilarityMapForModel(mol, get_fp, get_proba)
    imgdata = io.StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)  # rewind the data
    plt.savefig(imgdata, format="svg", bbox_inches="tight")

    return imgdata.getvalue()


def main(smiles, calculate_ad=False, make_prop_img=False, **kwargs):

    print(smiles)

    def default(key, d):
        if key in d.keys():
            return d[key]
        else:
            return False

    models = sorted([f for f in glob.glob(os.path.join(os.path.dirname(os.path.realpath(__file__)), "models/assay_inter/*.pkl"))])

    values = {}

    for model in models:
        if not default(MODEL_DICT[os.path.basename(model)], kwargs):
            continue
        with open(model, 'rb') as f:
            m = pickle.load(f)

        # for now I have no AD data so we are just doing this
        # pred, pred_proba, ad = run_prediction(m, None, smiles, calculate_ad=calculate_ad)
        pred, pred_proba, ad = run_prediction(m, None, smiles, calculate_ad=False)

        svg_str = ""
        if make_prop_img:
            svg_str = get_prob_map(m, smiles)

        values[MODEL_DICT[os.path.basename(model)]] = [OUTCOME_DICT[MODEL_DICT[os.path.basename(model)]][int(pred)], str(round(float(pred_proba) * 100, 2)) + "%", "NA", svg_str]

    return [[key] + val for key, val in values.items()]


def write_csv_file(smiles_list, calculate_ad=False):
    headers = list(MODEL_DICT.keys())

    if calculate_ad:
        headers = headers + [_ + "_AD" for _ in headers]

    string_file = StringIO()
    writer = csv.DictWriter(string_file, fieldnames=['SMILES', *headers])
    writer.writeheader()

    for smiles in tqdm(smiles_list):
        molecule = MolFromSmiles(smiles)

        row = {'SMILES': smiles}

        if molecule is None:
            row['SMILES'] = f"(invalid){smiles}"
            writer.writerow(row)
            continue

        data = main(smiles, calculate_ad=calculate_ad, **MODEL_DICT)

        for model_name, pred, pred_proba, ad, _ in data:
            try:
                pred_proba = float(pred_proba[:-1]) / 100  # covert back to 0-1 float
                row[
                    model_name] = pred_proba if pred == 1 else 1 - pred_proba  # this is to make sure its proba for class 1
            except ValueError:
                row[model_name] = "No prediction"  # if pred_proba is string skip
            if calculate_ad:
                row[model_name + "_AD"] = ad

        writer.writerow(row)

    return string_file.getvalue()


if __name__ == "__main__":
    import argparse
    import csv
    from io import StringIO
    from rdkit.Chem import MolFromSmiles
    from tqdm import tqdm

    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True,
                        help="location to csv of SMILES")
    parser.add_argument("--outfile", type=str, default=os.path.join(os.getcwd(), "phakin_output.csv"),
                        help="location and file name for output")
    parser.add_argument("--smiles_col", type=str, default="SMILES",
                        help="column name containing SMILES of interest"),
    parser.add_argument("--ad", action="store_true",
                        help="calculate the AD")
