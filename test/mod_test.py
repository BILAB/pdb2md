#%%
from modeller import *
from modeller.automodel import *

import Bio.PDB as pdb
import pymol2
import requests

#%%
def convert_complex_to_monomer(ID_dir="./",
                               pdb_path="./6m7f.pdb",
                               output_pdb_name="6m7f_mono"):
    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(pdb_path)
    save_dir = ID_dir + output_pdb_name +".pdb"
    pymol2_session.cmd.save(save_dir, "chain A and not resn HOH")
    pymol2_session.stop()

convert_complex_to_monomer()

#%%


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



#%%
def modelling_missing_res(id_dir="./",
                          pdb_id="6m7f",
                          pdb_path="./6m7f_mono.pdb"):

    print(f"Downloading {pdb_id}...")
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    data = requests.get(url).content
    if "not found" in str(data):
        print(f"Error: pdbID {pdb_id} is not found")
        pass
    else:
        with open(f"{id_dir}{pdb_id}.fasta", "wb") as f:
            f.write(data)
            f.close()

    # アラインメントをとる
    env = Environ()
    aln = Alignment(env)

    #mdlとfastaを読み込む順番でfastaよめなくなる？？？
    aln.append(file=f"6m7f.fasta",
               alignment_format="FASTA")

    mdl = Model(env, file=pdb_path)
    aln.append_model(mdl, align_codes=pdb_id)

    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,)
    aln.write(file="6m7f.ali", alignment_format="PIR")

    env = Environ()
    # HETATMの座標を保存するようにフラグを追加
    env.io.hetatm = True
    # 水分子の座標を保存するようにフラグを追加
    env.io.water = True

    a = AutoModel(env, alnfile='6m7f.ali', knowns='', sequence='6m7f')

    a.md_level = refine.fast

    a.starting_model = 1
    a.ending_model = 2
    a.make()

modelling_missing_res()

# %%
def download_full_fasta_for_pdb(pdb_id: str) -> str:
# %%
