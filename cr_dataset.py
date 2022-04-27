import argparse
from utils import *
import selfies as sf

def load_file(fname, augment):
    rxn_pairs = []
    with open(fname,"r") as f:
        for l in f.readlines():
            r,p=raw2pairs(l)
            if augment:
                aug = augment_rxn(r,p)
                rxn_pairs = rxn_pairs + aug 
            else:
                rxn_pairs.append( (r,p) )
    return rxn_pairs

def write_smiles(rxn_pairs, out_fname):
    with open(f"sm.{out_fname}.src", "w") as out_src, open(f"sm.{out_fname}.tgt", "w") as out_tgt:
        for r,p in rxn_pairs:
            out_src.write( smi_tokenizer( ".".join(r) ) + "\n")
            out_tgt.write( smi_tokenizer( ".".join(p) ) + "\n")

def write_deepsmiles(rxn_pairs, out_fname):
    with open(f"ds.{out_fname}.src", "w") as out_src, open(f"ds.{out_fname}.tgt", "w") as out_tgt:
        for r,p in rxn_pairs:
            r = map(sm2ds, r)
            p = map(sm2ds, p)
            out_src.write( smi_tokenizer( ".".join(r) ) + "\n")
            out_tgt.write( smi_tokenizer( ".".join(p) ) + "\n")

def write_selfies(rxn_pairs, out_fname):
    with open(f"sf.{out_fname}.src", "w") as out_src, open(f"sf.{out_fname}.tgt", "w") as out_tgt:
        for r,p in rxn_pairs:
            r = map(sf.encoder, r)
            p = map(sf.encoder, p)
            line_src = sf_tokenizer( ".".join(r) ) + "\n"
            line_tgt = sf_tokenizer( ".".join(p) ) + "\n"
            # selfies lib has bug in it where it fails to generate proper
            # SMILES if line ends with expl. Filter these out
            if (line_src[-5:] == "expl\n") or (line_tgt[-5:] == "expl\n"):
                continue
            out_src.write(line_src)
            out_tgt.write(line_tgt)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str)
    parser.add_argument("--output", type=str)
    parser.add_argument("--augment", type=int)
    
    args = parser.parse_args()

    rp = load_file(args.input, args.augment == 1)
    write_smiles(rp, args.output)
    write_deepsmiles(rp, args.output)
    write_selfies(rp, args.output)

main()

