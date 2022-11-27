# %%
import os
import shutil
import configparser
import subprocess
import Bio.PDB as pdb

def make_dirs_for_results(id_list: list) -> list:
    pwd = os.path.abspath(os.path.dirname(__file__))
    pdb_dir = os.path.join(pwd, "init_pdb")
    dihed_dir = os.path.join(pwd, "dihed")
    dists_dir = os.path.join(pwd, "dists")
    os.makedirs(name=pdb_dir, exist_ok=True)
    os.makedirs(name=dihed_dir, exist_ok=True)
    os.makedirs(name=dists_dir, exist_ok=True)
    return [pdb_dir, dihed_dir, dists_dir]

def make_init_files(id_list: list) -> None:
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
            if r.get_resname() == "FPP":
                fppid = str(r.get_id()[1])

        dihed_script = open(f"{id}/amber/pr/dihed.in", "w")
        dihed_script.write(f"""trajin ./001/mdcrd 1 last 50
trajin ./002/mdcrd 1 last 50
dihedral dih1 :{fppid}@O1 :{fppid}@C1 :{fppid}@C2 :{fppid}@C3 out dihed.txt range360
strip :SOD,CLA,WAT,TIP3,Na,Cl-
""")
        dihed_script.close()
        four_dist_script = open(f"{id}/amber/pr/4dist.in", "w")
        four_dist_script.write(f"""trajin ./001/mdcrd 1 last 50
trajin ./002/mdcrd 1 last 50
# C1toC6
distance dist1 :{fppid}@C1 :{fppid}@C7 out 4dist.txt
# C1toC7
distance dist2 :{fppid}@C1 :{fppid}@C8 out 4dist.txt
# C1toC10
distance dist3 :{fppid}@C1 :{fppid}@C12 out 4dist.txt
# C1toC11
distance dist4 :{fppid}@C1 :{fppid}@C13 out 4dist.txt
""")
        four_dist_script.close()

config = configparser.ConfigParser(strict=False, allow_no_value=True)
config.read("config.ini")
id = [k.upper() for k, v in config.items("ID")]
pdb_dir, dihed_dir, dists_dir = make_dirs_for_results(id_list=id)

for i in id:
    get_init_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./trajfix.in -p ../../top/leap.parm7
cp ./init.pdb {pdb_dir}/{i}.pdb
cp ./traj.trr {pdb_dir}/{i}.trr
cd ../../../
"""
    subprocess.run(get_init_cmd, shell=True)

make_init_files(id_list=id)

for i in id:
    get_dihed_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./dihed.in -p ../../top/leap.parm7
cp ./dihed.txt {dihed_dir}/{i}.txt
cd ../../../
"""
    subprocess.run(get_dihed_cmd, shell=True)

    get_dists_cmd = f"""
cd ./{i}/amber/pr
/home/apps/amber22/bin/cpptraj -i ./4dist.in -p ../../top/leap.parm7
cp ./4dist.txt {dists_dir}/{i}.txt
cd ../../../
"""
    subprocess.run(get_dists_cmd, shell=True)