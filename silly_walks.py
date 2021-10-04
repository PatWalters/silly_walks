# 1/usr/bin/env python

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd


class SillyWalks:
    def __init__(self, df):
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
            silly_bits = [
                x for x in [self.count_dict.get(x) for x in on_bits] if x is None
            ]
            score = len(silly_bits) / len(on_bits)
        else:
            score = 1
        return score


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Silly Walks")

    parser.add_argument("smiles", type=str, help="CSV file with SMILES")
    parser.add_argument(
        "-s", "--smiles_col", default="SMILES", type=str, help="SMILES column name"
    )
    parser.add_argument(
        "--sep", default=",", type=str, help="Column separator in input file"
    )
    parser.add_argument(
        "-r", "--reference", default="chembl_drugs.smi", type=str, help="Reference file"
    )
    parser.add_argument(
        "--strip_out",
        action="store_true",
        help="Output only SMILES and sillyness score",
    )
    parser.add_argument("-o", "--out", default="silly.csv", type=str, help="Output")

    args = parser.parse_args()

    df = pd.read_csv(args.smiles, sep=args.sep)
    df_ref = pd.read_csv(
        args.reference, sep=" ", names=["SMILES", "Name"]
    )  # TODO: Generalize?

    silly_walks = SillyWalks(df_ref)

    df["silly"] = df[args.smiles_col].apply(silly_walks.score)
    df.sort_values("silly", ascending=False, inplace=True)

    if args.strip_out:
        df[[args.smiles_col, "silly"]].to_csv(args.out, index=False)
    else:
        df.to_csv(args.out, index=False, float_format="%.3f")
