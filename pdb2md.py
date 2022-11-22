import os
import requests
import configparser
import shutil
import subprocess
import Bio.PDB as PDB
from Bio import SeqIO
from modeller import *
from modeller.automodel import *
import pymol2
import metapredict as meta
from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.protocols import minimization_packing as pack_min
from absl import flags

def path_to_abspath(path: str) -> str:
    if path:
        if "~" in path:
            path = os.path.expanduser(path)
        if "./" in path:
            path = path.replace("./",
                                os.path.abspath(".") + os.sep)
        path = os.path.normpath(path)
    return path

def make_id_dirs(ID_list: list,
                 dist_dir: str,
                 dir_name: str) -> dict:
    dist_dir = path_to_abspath(dist_dir)
    workbench_dir_path = f"{dist_dir}/{dir_name}"
    os.makedirs(workbench_dir_path,
                exist_ok=True)
    id_dirs = {}
    for ID in ID_list:
        id_dir = os.path.join(workbench_dir_path,
                              ID)
        os.makedirs(id_dir,
                    exist_ok=True)
        id_dirs[ID] = id_dir
    return id_dirs

def insert_templete_residue(residue_name: str,
                            output_pdb_name: str,
                            pdb_path: str,
                            templete_pdb_path: str) -> None:
    pdb_path = path_to_abspath(pdb_path)
    templete_pdb_path = path_to_abspath(templete_pdb_path)
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
                               output_pdb_name: str) -> None:
    pdb_path = path_to_abspath(pdb_path)
    PDBparser = PDB.PDBParser()
    seq = ""
    struct = PDBparser.get_structure("seq",
                                     pdb_path)
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
        shutil.copy(pdb_path,
                    "temp.pdb")
        os.rename("temp.pdb",
                  output_pdb_name)

def preparemd_settings_to_list(config: configparser.ConfigParser,
                               pdb_path: str,
                               distdir: str,
                               preparemd_path: str) -> list:
    preparemd_config = config.items("PREPAREMD_SETTINGS")
    distdir = path_to_abspath(distdir)
    pdb_path = path_to_abspath(pdb_path)
    preparemd_path = path_to_abspath(preparemd_path)
    mol2_name = config["RESIDUES_NAME_IN_TEMPLETE"]["insert_substrate_name"]
    mol2_path = config["PATH"]["parameter_file_path"]
    mol2_path = path_to_abspath(mol2_path)
    preparemd_cmd = ["python3",
                     preparemd_path,
                     "--file",
                     f"{pdb_path}",
                     "--distdir",
                     f"{distdir}",
                     "--mol2",
                     f"{mol2_name} = loadMol2 {mol2_path}"]
    for k, v in preparemd_config:
        preparemd_cmd.append(f"--{k}")
        preparemd_cmd.append(v)
    for i, v in enumerate(preparemd_cmd):
        preparemd_cmd[i] = path_to_abspath(preparemd_cmd[i])
    return preparemd_cmd

def rosetta_packing_residues(pdb_path: str,
                             output_pdb_dir: str,
                             output_pdb_name: str) -> None:
    pdb_path = path_to_abspath(pdb_path)
    output_pdb_dir = path_to_abspath(output_pdb_dir)

    init_options = ""
    for k, v in config.items("ROSETTA_SETTINGS"):
        init_options += f"-{k} {v} "
    init(init_options)

    pose = pose_from_pdb(pdb_path)

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())

    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)

    if not os.getenv("DEBUG"):
        packer.apply(pose)

    pose.dump_pdb(f"{output_pdb_dir}/{output_pdb_name}")

    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(f"{output_pdb_dir}/{output_pdb_name}", "packed")
    pymol2_session.cmd.select("H", "hydro")
    pymol2_session.cmd.save(F"{output_pdb_dir}/{output_pdb_name}", "packed and not H")
    pymol2_session.stop()

def make_qscript(workbench_dir: str,
                 id_dirs: dict) -> None:
    dir_list_for_sh_script = ""
    for ID, dir in id_dirs.items():
        dir_list_for_sh_script += f"\t{os.path.basename(dir)}\n"
    qsub_sh_path = f"{workbench_dir}/qsub.sh"
    qsub_sh_path = path_to_abspath(qsub_sh_path)

    with open(qsub_sh_path,"w") as f:
        f.write(f"""#!/bin/zsh
DIR=(\n{dir_list_for_sh_script})
for i in $DIR
do
    echo $1
    cd ./$i/amber
    /usr/local/bin/qsub ./totalrun.sh -N $i
    cd ../../
done""")
    f.close()

def download_pdb_files(pdb_id: str,
                       id_dir: str) -> str:
    id_dir = path_to_abspath(id_dir)
    if len(pdb_id) == 4:
        print(f"Downloading {pdb_id}...")
        url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
        data = requests.get(url).content
        if "not found" in str(data):
            print(f"Error: pdbID {pdb_id} is not found")
            pass
        with open(f"{id_dir}/{pdb_id}.pdb", "wb") as f:
            f.write(data)
            f.close()
            return os.path.abspath(f.name)
    else:
        print(f"Downloading {pdb_id}...")
        url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v3.pdb"
        data = requests.get(url).content
        if "not exist" in str(data):
            print(f"Error: uniplotID {pdb_id} is not found")
            pass
        with open(f"{id_dir}/AF-{pdb_id}-F1-model_v3.pdb", "wb") as f:
            f.write(data)
            f.close()
            return os.path.abspath(f.name)

def remove_alreadyexist_workbench(workbench_dir: str,
                                  flag: bool) -> None:
    if os.path.isdir(workbench_dir):
        print(f"{workbench_dir} already exists")
        if flag == True:
            shutil.rmtree(workbench_dir)
            print(f"Removed {workbench_dir}.")

def convert_complex_to_monomer(id_dir: str,
                               pdb_path: str,
                               output_pdb_name: str) -> str:
    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(pdb_path)
    pymol2_session.cmd.save(os.path.join(id_dir,
                                         output_pdb_name),
                            "chain A and not resn HOH",
                            format="pdb")
    pymol2_session.stop()
    return os.path.join(id_dir, output_pdb_name)

def download_fasta(pdb_id: str,
                   modeller_dir: str) -> str:
    modeller_dir = path_to_abspath(modeller_dir)
    print(f"Downloading FASTA file of {pdb_id}...")
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}"
    data = requests.get(url).content
    if "not found" in str(data):
        print(f"Error: pdbID {pdb_id} is not found")
        pass
    else:
        with open(f"{modeller_dir}/{pdb_id}.fasta", "wb") as f:
            f.write(data)
            f.close()
    fasta_path = f"{modeller_dir}/{pdb_id}.fasta"
    fasta_path = path_to_abspath(fasta_path)
    return fasta_path

def change_align_code_in_fasta(fasta_path: str,
                               align_codes: str,
                               alignment_format: str) -> str:
    fasta_path = path_to_abspath(fasta_path)
    for record in SeqIO.parse(fasta_path,
                              alignment_format):
        id = record.id
        desc = record.description
        seq = record.seq

    record.id = align_codes
    SeqIO.write(record,
                fasta_path,
                alignment_format)
    return align_codes

def modelling_missing_res(pdb_id: str,
                          id_dir: str,
                          modeller_dir: str,
                          fasta_path: str,
                          pdb_path: str) -> str:

    modeller_dir = path_to_abspath(modeller_dir)
    fasta_path = path_to_abspath(fasta_path)
    pdb_path = path_to_abspath(pdb_path)

    code_fill = change_align_code_in_fasta(fasta_path=fasta_path,
                                           align_codes=f"{pdb_id}_fill",
                                           alignment_format="fasta")

    env = Environ()
    aln = Alignment(env)
    env.io.atom_files_directory = [modeller_dir]

    aln.append(file=fasta_path,
               align_codes=code_fill,
               alignment_format="fasta")
    mdl = Model(env,
                file=pdb_path)
    aln.append_model(mdl,
                     align_codes=pdb_id)
    aln.salign(rms_cutoff=3.5,
               normalize_pp_scores=False)
    aln.write(file=f"{modeller_dir}/{pdb_id}.ali",
              alignment_format="pir")

    env = Environ()
    env.io.atom_files_directory = [id_dir]
    env.io.hetatm = True
    env.io.water = True
    a = AutoModel(env,
                  alnfile=f'{modeller_dir}/{pdb_id}.ali',
                  knowns=pdb_id,
                  sequence=code_fill)
    a.md_level = refine.fast
    a.starting_model = 1
    a.ending_model = 1
    a.make()

    shutil.move(f"./{code_fill}.B99990001.pdb",
                f"{modeller_dir}/{code_fill}.B99990001.pdb")
    shutil.move(f"./{code_fill}.D00000001",
                f"{modeller_dir}/{code_fill}.D00000001")
    shutil.move(f"./{code_fill}.ini",
                f"{modeller_dir}/{code_fill}.ini")
    shutil.move(f"./{code_fill}.rsr",
                f"{modeller_dir}/{code_fill}.rsr")
    shutil.move(f"./{code_fill}.sch",
                f"{modeller_dir}/{code_fill}.sch")
    shutil.move(f"./{code_fill}.V99990001",
                f"{modeller_dir}/{code_fill}.V99990001")

    pdb_path_modelled = f"{modeller_dir}/{code_fill}.B99990001.pdb"

    return pdb_path_modelled

config_path = "./config.ini"

flags.DEFINE_string(name="config_file",
                    default=config_path,
                    help="config file path. default is ./config.ini")
FLAGS = flags.FLAGS
if __name__ == "__main__":
    FLAGS(sys.argv)
    config_path = FLAGS.config_file

config = configparser.ConfigParser(allow_no_value=True,
                                   strict=False,
                                   delimiters="=")
config.optionxform = str
config.read(config_path)

flg_remove_disordered_residue: bool = config.getboolean("SETTINGS",
                                                        "remove_disordered_residue")
flg_insert_residue_from_templete: bool = config.getboolean("SETTINGS",
                                                           "insert_residue_from_temolete")
flg_rosetta_packing: bool = config.getboolean("SETTINGS",
                                              "rosetta_packing")
flg_insert_substrate_from_templete: bool = config.getboolean("SETTINGS",
                                                             "insert_substrate_from_templete")
flg_make_brandnew_workbench_if_existed: bool = config.getboolean("SETTINGS",
                                                                 "make_brandnew_workbench_if_existed")
flg_convert_complex_to_monomer: bool = config.getboolean("SETTINGS",
                                                         "convert_complex_to_monomer")
flg_modelling_missing_residue: bool = config.getboolean("SETTINGS",
                                                        "modelling_missing_residue")

workbench_dir = os.path.join(config["PATH"]["distination_path"],
                             config["SETTINGS"]["workbench_dir_name"])
workbench_dir = path_to_abspath(workbench_dir)

remove_alreadyexist_workbench(workbench_dir=workbench_dir,
                              flag=flg_make_brandnew_workbench_if_existed)

ID_list: list = [key.upper() for key in config["ID"]]
id_dirs: dict = make_id_dirs(ID_list=ID_list,
                             dist_dir=config["PATH"]["distination_path"],
                             dir_name=config["SETTINGS"]["workbench_dir_name"])
ID_pdb_paths: dict = {}

for ID, dir in id_dirs.items():
    ID_pdb_paths[ID] = download_pdb_files(pdb_id=ID,
                                          id_dir=dir)

if flg_convert_complex_to_monomer == True:
    for ID, pdb_path in ID_pdb_paths.items():
        if len(ID) == 4:
            print(f"Converting complex {ID} to monomer...")
            ID_pdb_paths[ID] = convert_complex_to_monomer(id_dir=id_dirs[ID],
                                                          pdb_path=pdb_path,
                                                          output_pdb_name=f"{ID}_monomer.pdb")

if flg_modelling_missing_residue == True:
    for ID, pdb_path in ID_pdb_paths.items():
        if len(ID) == 4:
            modeller_dir = os.path.join(id_dirs[ID],
                                        "modeller")
            os.mkdir(modeller_dir)
            fasta_path = download_fasta(pdb_id=ID,
                                        modeller_dir=modeller_dir)
            print(f"Modelling missing residue of {ID}...")
            ID_pdb_paths[ID] = modelling_missing_res(pdb_id=ID,
                                                     id_dir=id_dirs[ID],
                                                     modeller_dir=modeller_dir,
                                                     fasta_path=fasta_path,
                                                     pdb_path=pdb_path)

if flg_remove_disordered_residue == True:
    for ID, pdb_path in ID_pdb_paths.items():
        print(f"Removing disordered residues from {ID}...")
        remove_disordered_residues(pdb_path=pdb_path,
                                    output_pdb_name=f"{id_dirs[ID]}/{ID}_disordered_removed.pdb")
        ID_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_disordered_removed.pdb"
else:
    print(f"Skip removing disordered residues")

if flg_insert_residue_from_templete == True:
    for ID, pdb_path in ID_pdb_paths.items():
        templete_residue_name = config["RESIDUES_NAME_IN_TEMPLETE"]["insert_residue_name"]
        print(f"Inserting {templete_residue_name} into {ID} from templete...")
        insert_templete_residue(residue_name=templete_residue_name,
                                output_pdb_name=f"{id_dirs[ID]}/{ID}_res_inserted.pdb",
                                pdb_path=pdb_path,
                                templete_pdb_path=config["PATH"]["templete_pdb_path"])
        ID_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_res_inserted.pdb"
else:
    print(f"Inserting residues from templete is skipped")

if flg_rosetta_packing == True:
    for ID, pdb_path in ID_pdb_paths.items():
        print(f"Packing {ID} for optimizing sidechains and hetero metals...")
        rosetta_packing_residues(pdb_path=pdb_path,
                                output_pdb_dir=id_dirs[ID],
                                output_pdb_name=f"{ID}_packed.pdb")
        ID_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_packed.pdb"
else:
    print(f"Skip packing...")

if flg_insert_substrate_from_templete == True:
    substrate_name = config["RESIDUES_NAME_IN_TEMPLETE"]["insert_substrate_name"]
    for ID, pdb_path in ID_pdb_paths.items():
        print(f"Inserting {substrate_name} into {ID} from templete...")
        insert_templete_residue(residue_name=substrate_name,
                                output_pdb_name=f"{id_dirs[ID]}/{ID}_res_sub_inserted.pdb",
                                pdb_path=pdb_path,
                                templete_pdb_path=config["PATH"]["templete_pdb_path"])
        ID_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_res_sub_inserted.pdb"
else:
    print(f"Inserting substrate from templete is skipped")

for ID, pdb_path in ID_pdb_paths.items():
    print(f"Excuting preparemd.py for {ID}...")
    preparemd_cmd = preparemd_settings_to_list(config=config,
                                               pdb_path=pdb_path,
                                               distdir=id_dirs[ID],
                                               preparemd_path=config["PATH"]["preparemd_script_path"])
    subprocess.run(preparemd_cmd)

shutil.copy("config.ini", workbench_dir)
print("Copied config.ini to workbench directory.")

print("Making a script file for submitting MD jobs in yayoi...")
make_qscript(workbench_dir=config["PATH"]["distination_path"] +
                           "/" +
                           config["SETTINGS"]["workbench_dir_name"],
             id_dirs=id_dirs)

print("Process terminated.")
