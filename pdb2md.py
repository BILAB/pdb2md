# %%
import configparser
import os
import requests
import Bio.PDB as PDB
import pymol2
import metapredict as meta
import shutil

config = configparser.ConfigParser(allow_no_value=True, strict=False)
config.read("config.ini")

def make_ID_dirs(ID_list: list,
                 name: str):
    parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), name))
    os.makedirs(parent_dir, exist_ok=True)
    ID_dirs = {}
    for ID in ID_list:
        ID_dir = os.path.join(parent_dir, ID)
        os.makedirs(ID_dir, exist_ok=True)
        ID_dirs[ID] = ID_dir
    return ID_dirs

def insert_templete_residue(residue_name: str,
                            output_pdb_name: str,
                            pdb_path: str,
                            templete_pdb_path: str = config["PATH"]["templete_pdb_path"]):
    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(templete_pdb_path, "templete")
    pymol2_session.cmd.load(pdb_path, "objective")
    pymol2_session.cmd.super("templete", "objective")
    pymol2_session.cmd.select("res", f"templete and resn {residue_name}")
    pymol2_session.cmd.extract("res", "res")
    pymol2_session.cmd.delete("templete")
    pymol2_session.cmd.save(output_pdb_name, "res or objective")
    pymol2_session.stop()

def res3to1(res: str) -> str:
    if res == "ALA":
        return "A"
    elif res == "ARG":
        return "R"
    elif res == "ASN":
        return "N"
    elif res == "ASP":
        return "D"
    elif res == "CYS":
        return "C"
    elif res == "GLN":
        return "Q"
    elif res == "GLU":
        return "E"
    elif res == "GLY":
        return "G"
    elif res == "HIS":
        return "H"
    elif res == "ILE":
        return "I"
    elif res == "LEU":
        return "L"
    elif res == "LYS":
        return "K"
    elif res == "MET":
        return "M"
    elif res == "PHE":
        return "F"
    elif res == "PRO":
        return "P"
    elif res == "SER":
        return "S"
    elif res == "THR":
        return "T"
    elif res == "TRP":
        return "W"
    elif res == "TYR":
        return "Y"
    elif res == "VAL":
        return "V"
    else:
        return None

def remove_disordered_residues(pdb_path: str,
                               output_pdb_name: str):
    PDBparser = PDB.PDBParser()
    seq = ""
    struct = PDBparser.get_structure("seq", pdb_path)
    res = struct[0]["A"].get_list()
    for i in res:
        resname = res3to1(i.resname)
        if resname != None:
            seq += resname
        else:
            print(f"Unknown residue: {i.id[1]} {i.resname} in {pdb_path}.")

    metapredict = meta.predict_disorder(seq)
    metapredict_domains = meta.predict_disorder_domains(seq).disordered_domain_boundaries

    remove_cmd = ""
    # metapredict_domainsに要素がある場合
    if len(metapredict_domains) != 0:
        for metapredict_domain in metapredict_domains:
            remove_cmd += f"resi {metapredict_domain[0]}:{metapredict_domain[1]} or "
        remove_cmd = remove_cmd[:-4]

        pymol2_session = pymol2.PyMOL()
        pymol2_session.start()
        pymol2_session.cmd.load(pdb_path, "master")
        pymol2_session.cmd.select("disordered", remove_cmd)
        pymol2_session.cmd.save(output_pdb_name, "master and not disordered")
        pymol2_session.stop()
    else:
        shutil.copy(pdb_path, "temp.pdb")
        os.rename("temp.pdb", output_pdb_name)

ID_list = [key.upper() for key in config["ID"]]
ID_dirs = make_ID_dirs(ID_list, "pdb2md")
ID_pdb_paths = {}

for ID, dir in ID_dirs.items():
    #IDが4文字の場合
    if len(ID) == 4:
        url = "https://files.rcsb.org/download/" + str.upper(ID) + ".pdb"
        data = requests.get(url).content
        if "not exist" in str(data):
            print(f"Error: pdbID {(ID)} is not found")
            continue
        with open(f"{dir}/{ID}.pdb", "wb") as f:
            f.write(data)
            ID_pdb_paths[ID] = os.path.abspath(f.name)
    else:
        url = "https://alphafold.ebi.ac.uk/files/AF-" + ID + "-F1-model_v3.pdb"
        data = requests.get(url).content
        if "not exist" in str(data):
            print(f"Error: uniplotID {ID} is not found")
            continue
        with open(f"{dir}/AF-{ID}-F1-model_v3.pdb", "wb") as f:
            f.write(data)
            ID_pdb_paths[ID] = os.path.abspath(f.name)

for ID, pdb_path in ID_pdb_paths.items():
    if len(ID) != 4:
        remove_disordered_residues(pdb_path, f"{ID_dirs[ID]}/{ID}_disordered_removed.pdb")
