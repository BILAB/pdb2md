#%%
from modeller import *
from modeller.automodel import *
import requests
from Bio import SeqIO

#%%
def download_fasta(pdb_id,
                   id_dir):
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
    return f

def change_align_code_in_fasta(file,
                               align_codes,
                               alignment_format):
    # fastaファイルのalign_codesを変更
    for record in SeqIO.parse("6m7f.fasta", 'fasta'):
        id = record.id
        desc = record.description
        seq = record.seq

    # idを変更
    record.id = "6m7f_fill"
    # 書き込み
    SeqIO.write(record, "6m7f.fasta", "fasta")

    return align_codes

def modelling_missing_res(id_dir="./",
                          pdb_id="6m7f",
                          pdb_path="./6m7f_mono.pdb"):

    f = download_fasta(pdb_id=pdb_id,
                       id_dir=id_dir)

    env = Environ()
    aln = Alignment(env)
    code_fill = change_align_code_in_fasta(file="6m7f.fasta",
                                           align_codes="6m7f_fill",
                                           alignment_format="fasta")

    # mdlとfastaを読み込む順番でfastaよめなくなる？？？
    # fastaファイルをalnに読み込む
    aln.append(file="6m7f.fasta",
               align_codes="all",
               alignment_format="FASTA")
    # pdbファイルをalnに読み込む
    mdl = Model(env,
                file=pdb_path)
    aln.append_model(mdl,
                     align_codes=pdb_id)
    # アラインメントを取る
    aln.salign(rms_cutoff=3.5,
               normalize_pp_scores=False)
    # aliファイルに出力
    aln.write(file="6m7f.ali",
              alignment_format="PIR")

    # model構築のために初期化
    env = Environ()
    env.io.hetatm = True  # HETATMの座標を保存するようにフラグを追加
    env.io.water = True  # 水分子の座標を保存するようにフラグを追加

    # モデル構築の実行
    a = AutoModel(env,
                  alnfile='6m7f.ali',
                  knowns='6m7f',
                  sequence=code_fill)
    a.md_level = refine.fast
    a.starting_model = 1
    a.ending_model = 1
    a.make()

modelling_missing_res()

#%%
# 6m7f.fastaのrecord.idを書き換える


# %%
