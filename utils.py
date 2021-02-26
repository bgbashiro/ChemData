import deepsmiles as ds
import selfies as sf
from rdkit import Chem
import re

def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    pattern =  r"(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    # assert smi == ''.join(tokens)
    return ' '.join(tokens)

def sf_tokenizer(sfi):
    return ' '.join("." if tok=="." else tok[1:-1] for tok in sf.split_selfies(sfi))

def strip_atom_numbering(mol_str):
    mol = Chem.MolFromSmiles(mol_str)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)

def raw2pairs(sm_line):
    """
    Convert raw line from NIPS into (reagents, products) splits.
    Each reagent, products are list of molecules. The atom numbering is
    stripped
    """
    reagents, products = sm_line.split(" ")[0].split(">>")
    reagents = [ strip_atom_numbering(mol_str) for mol_str in  reagents.split(".")]
    products = [ strip_atom_numbering(mol_str) for mol_str in  products.split(".")]
    return reagents, products

def sm2ds(line):
    # Takes schwaller's preprocessed SMILES and turns them into deepSMILES
    converter = ds.Converter(rings=True, branches=True)
    line = line.replace(" ", "")
    molecules = line.split(".")
    new_line = []
    for molecule in molecules:
        new_molecule = converter.encode(molecule)
        new_line.append(new_molecule)
    new_line = ".".join(new_line)
    return new_line


def ds2sm(line):
    # Takes a deepSMILES line and turns it into SMILES
    converter = ds.Converter(rings=True, branches=True)
    line = line.replace(" ", "")
    molecules = line.split(".")
    new_line = []
    for molecule in molecules:
        new_molecule = converter.decode(molecule)
        new_line.append(new_molecule)
    new_line = ".".join(new_line)
    return new_line

def sm2sf(line):
    pass

def sf2sm(line):
    words = line.split(".")
    words = [ "[" + word.replace(" ", "][") + "]" for word in words ] 
    new_line = []
    for word in words: 
        new_line.append( sf.decoder(word))
    new_line = ".".join(new_line)
    return new_line

