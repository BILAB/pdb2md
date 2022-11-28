# %%
import os
import shutil
import configparser
import subprocess
import Bio.PDB as pdb


def make_dirs_for_results(id_list: list) -> list:
    """測定結果を格納するディレクトリを作成する

    Args:
        id_list (list): config.iniなどから読み込んだIDのリスト

    Returns:
        list: 指定したIDの解析結果を格納する用のディレクトリリスト
    """

    pwd = os.path.abspath(os.path.dirname(__file__))
    pdb_dir = os.path.join(pwd, "init_pdb")
    dihed_dir = os.path.join(pwd, "dihed")
    dists_dir = os.path.join(pwd, "dists")
    os.makedirs(name=pdb_dir, exist_ok=True)
    os.makedirs(name=dihed_dir, exist_ok=True)
    os.makedirs(name=dists_dir, exist_ok=True)
    return [pdb_dir, dihed_dir, dists_dir]


def dihed_to_in_file(dihed_dict: dict,
                     resi: str,
                     id: str):
    """辞書を読み取り、二面角測定用のinputファイルを作成する

    Args:
        dists_dict (dict): 測定したい距離の辞書 (key: 結果名, value: 距離を構成する2原子のリスト)
        resi (str): 二面角を測定したい残基
        id (str): 対象となるタンパク質のID
    """

    dihed_script = open(f"{id}/amber/pr/dihed.in", "w")
    dihed_script.write(f"""trajin ./001/mdcrd 1 last 50
trajin ./002/mdcrd 1 last 50
""")
    dihed_script.close()
    for name, atoms in dihed_dict.items():
        dihed_script = open(f"{id}/amber/pr/dihed.in", "a")
        dihed_script.write(f"""dihedral {name} :{resi}@{atoms[0]} :{resi}@{atoms[1]} :{resi}@{atoms[2]} :{resi}@{atoms[3]} out dihed.txt range360
strip :SOD,CLA,WAT,TIP3,Na,Cl-
""")
        dihed_script.close()


def dists_to_in_file(dists_dict: dict,
                     resi: str,
                     id: str):
    """辞書を読み取り、距離測定用のinputファイルを作成する

    Args:
        dists_dict (dict): 測定したい距離の辞書 (key: 結果名, value: 距離を構成する2原子のリスト)
        resi (str): 距離を測定したい残基
        id (str): 対象となるタンパク質のID
    """

    dists_script = open(f"{id}/amber/pr/dists.in", "w")
    dists_script.write(f"""trajin ./001/mdcrd 1 last 50
trajin ./002/mdcrd 1 last 50
""")
    dists_script.close()
    for name, atoms in dists_dict.items():
        dists_script = open(f"{id}/amber/pr/dists.in", "a")
        dists_script.write(f"""distance {name} :{resi}@{atoms[0]} :{resi}@{atoms[1]} out dists.txt
""")
        dists_script.close()
    dists_script.close()


def make_in_files(id_list: list,
                  dihed_dict: dict,
                  dists_dict: dict,
                  res_name: str) -> None:
    """初期構造pdbファイルから対象残基の番号を取得し、二面角と距離測定用のinputファイルを作成する

    Args:
        id_list (list): config.iniなどから読み込んだIDのリスト
        dists_dict (dict): 測定したい距離の辞書 (key: 結果名, value: 距離を構成する2原子のリスト)
        dists_dict (dict): 測定したい距離の辞書 (key: 結果名, value: 距離を構成する2原子のリスト)
        res_name (str): 測定したい残基名
    """

    for id in id_list:
        pdbfile = os.path.join("./", id, "amber", "pr", "init.pdb")
        pdb_parser = pdb.PDBParser()
        try:
            struct = pdb_parser.get_structure("all", pdbfile)
        except FileNotFoundError:
            continue
        chain = struct[0].get_list()[0]
        res = chain.get_list()
        for r in res:
            if r.get_resname() == res_name:
                resi = str(r.get_id()[1])

        dihed_to_in_file(dihed_dict=dihed_dict,
                         resi=resi,
                         id=id)

        dists_to_in_file(dists_dict=dists_dict,
                         resi=resi,
                         id=id)


# config.iniの読み込み
config = configparser.ConfigParser(strict=False, allow_no_value=True)
config.read("config.ini")
id = [k.upper() for k, v in config.items("ID")]
pdb_dir, dihed_dir, dists_dir = make_dirs_for_results(id_list=id)

# pdbの初期構造とトラジェクトリを取得
for i in id:
    get_init_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./trajfix.in -p ../../top/leap.parm7
cp ./init.pdb {pdb_dir}/{i}.pdb
cp ./traj.trr {pdb_dir}/{i}.trr
cd ../../../"""
    subprocess.run(get_init_cmd, shell=True)

# 解析対象を定義
dihed_dict = {"O1~C3": ["O1", "C1", "C2", "C3"],
              "C1~C4": ["C1", "C2", "C3", "C4"]}
dists_dict = {"C1-C7": ["C1", "C7"],
              "C1-C8": ["C1", "C8"],
              "C1-C12": ["C1", "C12"],
              "C1-C13": ["C1", "C13"]}
make_in_files(id_list=id,
              res_name="FPP",
              dihed_dict=dihed_dict,
              dists_dict=dists_dict)

# 解析の実行
for i in id:
    get_dihed_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./dihed.in -p ../../top/leap.parm7
cp ./dihed.txt {dihed_dir}/{i}.txt
cd ../../../"""
    subprocess.run(get_dihed_cmd, shell=True)

    get_dists_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./dists.in -p ../../top/leap.parm7
cp ./dists.txt {dists_dir}/{i}.txt
cd ../../../"""
    subprocess.run(get_dists_cmd, shell=True)