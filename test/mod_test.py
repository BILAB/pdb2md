#%%
from modeller import *
from modeller.automodel import *
import requests

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

    print(f"Downloading FASTA file of {pdb_id}...")
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    data = requests.get(url).content
    if "not found" in str(data):
        print(f"Error: pdbID {pdb_id} is not found")
        pass
    else:
        with open(f"{id_dir}{pdb_id}.fasta", "wb") as f:
            f.write(data)
            f.close()

    env = Environ()
    aln = Alignment(env)

    # mdlとfastaを読み込む順番でfastaよめなくなる？？？
    # fastaファイルをalnに読み込む
    aln.append(file=f"6m7f.fasta",
               alignment_format="FASTA",
               align_codes="all")
    # pdbファイルをalnに読み込む
    mdl = Model(env,
                file=pdb_path)
    aln.append_model(mdl,
                     align_codes=pdb_id)
    # アラインメントを取る
    aln.salign(rms_cutoff=3.5,
               normalize_pp_scores=False,)
    # aliファイルに出力
    aln.write(file="6m7f.ali",
              alignment_format="PIR")

    # model構築のために初期化
    env = Environ()
    env.io.hetatm = True  # HETATMの座標を保存するようにフラグを追加
    env.io.water = True  # 水分子の座標を保存するようにフラグを追加

    a = AutoModel(env,
                  alnfile='6m7f.ali',
                  knowns='6m7f',
                  sequence='6m7f_fill')
    a.md_level = refine.fast
    a.starting_model = 1
    a.ending_model = 1
    a.make()

modelling_missing_res()
