#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, unicode_literals
import argparse
from rdkit import Chem
import pandas as pd
from utils import ds2sm,sf2sm

def canonicalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=True) if mol is not None else ''

def get_rank(row, base, max_rank):
    return next(
        (
            i
            for i in range(1, max_rank + 1)
            if row['target'] == row[f'{base}{i}']
        ),
        0,
    )

def main(opt):
    with open(opt.targets, 'r') as f:
        if opt.mol_format == "smiles":
            targets = [''.join(line.strip().split(' ')) for line in f.readlines()]
        elif opt.mol_format == "deepsmiles":
            targets = [canonicalize_smiles(ds2sm(''.join(line.strip().split(' ')))) for line in f.readlines()]
        elif opt.mol_format == "selfies":
            targets = [canonicalize_smiles(sf2sm(line)) for line in f.readlines()]

    predictions = [[] for _ in range(opt.beam_size)]

    test_df = pd.DataFrame(targets)
    print(test_df.shape)
    test_df.columns = ['target']
    total = len(test_df)

    with open(opt.predictions, 'r') as f:
        if opt.mol_format == "smiles":
            for i, line in enumerate(f.readlines()):
                pred_smile = ''.join(line.strip().split(' '))
                predictions[i % opt.beam_size].append(pred_smile)
        elif opt.mol_format == "deepsmiles":
            for i, line in enumerate(f.readlines()):
                pred_smile = ds2sm(''.join(line.strip().split(' ')))
                predictions[i % opt.beam_size].append(pred_smile)
        elif opt.mol_format == "selfies":
            for i, line in enumerate(f.readlines()):
                pred_smile = sf2sm(line)
                predictions[i % opt.beam_size].append(pred_smile)

    print(len(predictions))
    for i, preds in enumerate(predictions):
        print(len(preds))
        test_df[f'prediction_{i + 1}'] = preds
        test_df[f'canonical_prediction_{i + 1}'] = test_df[
            f'prediction_{i + 1}'
        ].apply(lambda x: canonicalize_smiles(x))


    test_df['rank'] = test_df.apply(lambda row: get_rank(row, 'canonical_prediction_', opt.beam_size), axis=1)

    correct = 0

    for i in range(1, opt.beam_size+1):
        correct += (test_df['rank'] == i).sum()
        invalid_smiles = (test_df[f'canonical_prediction_{i}'] == '').sum()
        if opt.invalid_smiles:
            print('Top-{}: {:.1f}% || Invalid SMILES {:.2f}%'.format(i, correct/total*100,
                                                                     invalid_smiles/total*100))
        else:
            print('Top-{}: {:.1f}%'.format(i, correct / total * 100))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='score_predictions.py',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-beam_size', type=int, default=5,
                       help='Beam size')
    parser.add_argument('-invalid_smiles', action="store_true",
                       help='Show % of invalid SMILES')
    parser.add_argument('-predictions', type=str, default="",
                       help="Path to file containing the predictions")
    parser.add_argument('-targets', type=str, default="",
                       help="Path to file containing targets")
    parser.add_argument('-mol_format', type=str, default="smiles",
                        help="smiles/deepsmiles or selfies")

    opt = parser.parse_args()
    main(opt)
