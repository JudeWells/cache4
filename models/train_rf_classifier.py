"""
Created by Jude Wells 2023-04-13
Trains RF balanced RF classifier on cache4 hits for CBLB
and decoys taken from Enamine HLL.
Does scaffold train / test split
Results for n_estimators 10 (pure binary prediction) thresholds 0.5, 0.5:
MCC:  0.66
bAcc: 0.73
"""
import warnings

import numpy as np
from imblearn.ensemble import BalancedRandomForestClassifier
from nonconformist.acp import AggregatedCp
from nonconformist.cp import IcpClassifier
from nonconformist.nc import NcFactory
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.metrics import matthews_corrcoef, balanced_accuracy_score
from sklearn.neighbors import KNeighborsRegressor

warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

from utils.scaffold_split import load_train_test_data

def prepare_data_for_classifier(train, test, RADIUS=3, NBITS=2048):
    x_train = np.array([np.array(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS)) for smi in train["smiles"]])
    x_test = np.array([np.array(AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), RADIUS, nBits=NBITS)) for smi in test["smiles"]])
    y_train = train["target"].values
    y_test = test["target"].values
    return x_train, x_test, y_train, y_test

def evaluate_preds(preds, true, hit_threshold=0.5, decoy_threshold=0.5):
    binary_pred = np.logical_and(preds[:, 0] < decoy_threshold, preds[:, 1] > hit_threshold)
    print("MCC: ", matthews_corrcoef(true, binary_pred))
    print("bAcc:", balanced_accuracy_score(true, binary_pred))


if __name__=="__main__":
    train, test = load_train_test_data()
    x_train, x_test, y_train, y_test = prepare_data_for_classifier(train, test)
    model = BalancedRandomForestClassifier(n_estimators=10, random_state=42)
    nc = NcFactory.create_nc(model, normalizer_model=KNeighborsRegressor(n_neighbors=11))
    icp = IcpClassifier(nc)
    acp = AggregatedCp(icp)
    acp.fit(x_train, y_train)
    test_pred = acp.predict(x_test)
    evaluate_preds(test_pred, y_test)


