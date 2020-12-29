# 1/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


class SillyWalks:
    def __init__(self, ref_filename):
        df = pd.read_csv(ref_filename, sep=" ", names=["SMILES", "Name"])
        self.count_dict = {}
        for smi in df.SMILES:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                fp = AllChem.GetMorganFingerprint(mol, 2)
                for k, v in fp.GetNonzeroElements().items():
                    self.count_dict[k] = self.count_dict.get(k, 0) + v

    def score(self, smiles_in):
        mol = Chem.MolFromSmiles(smiles_in)
        if mol:
            fp = AllChem.GetMorganFingerprint(mol, 2)
            on_bits = fp.GetNonzeroElements().keys()
            silly_bits = [x for x in [self.count_dict.get(x) for x in on_bits] if x is None]
            score = len(silly_bits) / len(on_bits)
        else:
            score = 1
        return score


if __name__ == "__main__":
    infile_name = sys.argv[1]
    df = pd.read_csv(infile_name)
    silly_walks = SillyWalks("chembl_drugs.smi")
    df['silly'] = df.SMILES.apply(silly_walks.score)
    df.sort_values('silly', ascending=False, inplace=True)
    print(df.head())
