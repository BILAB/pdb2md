import os
import sys
import requests
import configparser
import shutil
import subprocess
from absl import flags
import Bio.PDB as PDB
from Bio import SeqIO
import pymol2
import metapredict as meta
import modeller
import modeller.automodel as automodel
import pyrosetta.rosetta.utility as util
import pyrosetta.io as io
from pyrosetta.rosetta.core.init import init
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols import minimization_packing as pack_min


def path_to_abspath(path: str) -> str:
    """Normalize directory names.

    Args:
        path (str): Path to be normalized

    Returns:
        str: Normalized path
    """
    if path:
        if "~" in path:
            path = os.path.expanduser(path)
        if "./" in path:
            path = path.replace("./", os.path.abspath(".") + os.sep)
        path = os.path.normpath(path)
    return path


def make_id_dirs(id_list: list, dist_dir: str, dir_name: str) -> dict:
    """Create directories in workbench with the name of the PDB or uniprot ID.

    Args:
        id_list (list): PDB or uniprot ID list
        dist_dir (str): Directory path to create workbench
        dir_name (str): Workbench name

    Returns:
        dict: Dictionary of ID and directory path
    """
    dist_dir = path_to_abspath(dist_dir)
    workbench_dir_path = f"{dist_dir}/{dir_name}"
    os.makedirs(workbench_dir_path, exist_ok=True)
    id_dirs = {}
    for ID in id_list:
        id_dir = os.path.join(workbench_dir_path, ID)
        os.makedirs(id_dir, exist_ok=True)
        id_dirs[ID] = id_dir
    return id_dirs


def insert_templete_residue(pdb_path: str, templete_pdb_path: str,
                            residue_name: str, output_pdb_name: str) -> str:
    """Copy a residue with a given name from a template protein based on structural alignment.

    Args:
        pdb_path (str): Path of the pdb file into which the residues are to be inserted
        templete_pdb_path (str): Path of the pdb file to be used as a template
        residue_name (str): Residue name in the template structure to be copied
        output_pdb_name (str): Name of the pdb file to be inserted residues

    Returns:
        str: Path of the pdb file to be used as a template
    """
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
    return output_pdb_name


def res3to1(res: str) -> str:
    """Convert amino acids in 3-letter notation to 1-letter notation

    Args:
        res (str): Amino acids in 3-letter notation

    Returns:
        str: Amino acids in 1-letter notation
    """
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


def remove_disordered_residues(pdb_path: str, output_pdb_name: str) -> str:
    """Remove disordered area using the metapredict package to reduce computation

    Args:
        pdb_path (str): Path of the pdb file to be removed disordered area
        output_pdb_name (str): Pdb file name after removing the disordered area

    Returns:
        str: Path of pdb file after removing the disordered area
    """
    pdb_path = path_to_abspath(pdb_path)
    PDBparser = PDB.PDBParser()
    seq = ""
    struct = PDBparser.get_structure("seq", pdb_path)
    res = struct[0]["A"].get_list()
    for i in res:
        resname = res3to1(i.resname)
        if resname is not None:
            seq += resname
        else:
            print(f"Unknown residue: {i.id[1]} {i.resname} in {pdb_path}.")
    metapredict_domains = meta.predict_disorder_domains(
        seq
    ).disordered_domain_boundaries
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
        shutil.copy(pdb_path, "temp.pdb")
        os.rename("temp.pdb", output_pdb_name)
    return output_pdb_name


def preparemd_settings_to_list(
    config: configparser.ConfigParser, pdb_path: str, distdir: str, preparemd_path: str
) -> list:
    """Parse the preparemd settings in config.ini

    Args:
        config (configparser.ConfigParser): ConfigParser object reading config.ini
        pdb_path (str): Path of the pdb file to preparemd
        distdir (str): Path of directory to output preparemd results
        preparemd_path (str): Path to preparemd scripts

    Returns:
        list: List of formatted preparemd arguments and their values
    """
    preparemd_config = config.items("PREPAREMD_SETTINGS")
    distdir = path_to_abspath(distdir)
    pdb_path = path_to_abspath(pdb_path)
    preparemd_path = path_to_abspath(preparemd_path)
    mol2_name = config["RESIDUES_NAME_IN_TEMPLETE"]["insert_substrate_name"]
    mol2_path = config["PATH"]["parameter_file_path"]
    mol2_path = path_to_abspath(mol2_path)
    preparemd_cmd = [
        "python3",
        preparemd_path,
        "--file",
        f"{pdb_path}",
        "--distdir",
        f"{distdir}",
        "--mol2",
        f"{mol2_name} = loadMol2 {mol2_path}",
    ]
    for k, v in preparemd_config:
        preparemd_cmd.append(f"--{k}")
        preparemd_cmd.append(v)
    for i, v in enumerate(preparemd_cmd):
        preparemd_cmd[i] = path_to_abspath(preparemd_cmd[i])
    return preparemd_cmd


def remove_hydrogen(pdb_path: str, output_pdb_dir: str, output_pdb_name: str) -> str:
    """Removing water from pdb files output from rosetta

    Args:
        pdb_path (str): Pdb file to remove water (HOH)
        output_pdb_dir (str): Path of directory where PDB files will be saved after water removal
        output_pdb_name (str): Name of PDB file after water removal

    Returns:
        str: Path of PDB file after removing water
    """
    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(pdb_path, "packed")
    pymol2_session.cmd.select("H", "hydro")
    pymol2_session.cmd.save(f"{output_pdb_dir}/{output_pdb_name}", "packed and not H")
    pymol2_session.stop()
    output_pdb_path = os.path.join(output_pdb_dir, output_pdb_name)
    return output_pdb_path


def rosetta_packing_residues(
    pdb_path: str, output_pdb_dir: str, output_pdb_name: str
) -> str:
    """Perform packing (residue side-chain structure optimization) with PyRosetta

    Args:
        pdb_path (str): Pdb file to be packed
        output_pdb_dir (str): Path of directory where PDB files will be saved after Rosetta packing
        output_pdb_name (str): Name of PDB file after Rosetta packing

    Returns:
        str: Path of PDB file after Rosetta packing
    """

    pdb_path = path_to_abspath(pdb_path)
    output_pdb_dir = path_to_abspath(output_pdb_dir)

    init_options = ""
    for k, v in config.items("ROSETTA_SETTINGS"):
        init_options += f"-{k}={v} "
    init_options = util.string_split_simple(init_options, " ")
    init(init_options)

    pose = io.pose_from_pdb(pdb_path)

    tf = TaskFactory()
    tf.push_back(operation.RestrictToRepacking())

    packer = pack_min.PackRotamersMover()
    packer.task_factory(tf)

    if not os.getenv("DEBUG"):
        packer.apply(pose)

    pose.dump_pdb(f"{output_pdb_dir}/{output_pdb_name}")

    pdb_path = remove_hydrogen(
        pdb_path=pdb_path, output_pdb_dir=output_pdb_dir, output_pdb_name=output_pdb_name
    )
    return pdb_path


def make_qscript(workbench_dir: str, id_dirs: dict) -> str:
    """Create a shell script (qsub.sh) to submit output by preparemd to yayoi at once

    Args:
        workbench_dir (str): Path of workbench directory
        id_dirs (dict): Dictionary of ID and directory path

    Returns:
        str: Path to output shell scripts
    """
    dir_list_for_sh_script = ""
    for ID, dir in id_dirs.items():
        dir_list_for_sh_script += f"\t{os.path.basename(dir)}\n"
    qsub_sh_path = f"{workbench_dir}/qsub.sh"
    qsub_sh_path = path_to_abspath(qsub_sh_path)

    with open(qsub_sh_path, "w") as f:
        f.write(
            f"""#!/bin/zsh
DIR=(\n{dir_list_for_sh_script})
for i in $DIR
do
    echo $1
    cd ./$i/amber
    /usr/local/bin/qsub ./totalrun.sh -N $i
    cd ../../
done"""
        )
    f.close()
    return qsub_sh_path


def make_initscript(workbench_dir: str, id_dirs: dict) -> str:
    """Create a shell script (init.sh) that outputs a batch of PDB files with initial structure from MD results

    Args:
        workbench_dir (str): Path of workbench directory
        id_dirs (dict): Dictionary of ID and directory path

    Returns:
        str: Path to output shell scripts
    """
    dir_list_for_sh_script = ""
    for ID, dir in id_dirs.items():
        dir_list_for_sh_script += f"\t{os.path.basename(dir)}\n"
    init_sh_path = f"{workbench_dir}/init.sh"
    init_sh_path = path_to_abspath(init_sh_path)

    with open(init_sh_path, "w") as f:
        f.write(
            f"""#!/bin/zsh
modele load amber22
DIR=(\n{dir_list_for_sh_script})
for i in $DIR
do
    echo $1
    cd ./$i/amber/pr
    cpptraj -i ./trajfix.in -p ../../top/leap.parm7
    cd ../../../
done"""
        )
    f.close()
    return init_sh_path


def download_pdb_files(pdb_id: str, id_dir: str) -> str:
    """Download PDB files of PDB (in case of 4-letter ID) or uniprot ID

    Args:
        pdb_id (str): PDB or uniprot ID
        id_dir (str): Dictionary of ID and directory path

    Returns:
        str: Path to pdb file downloaded
    """
    id_dir = path_to_abspath(id_dir)
    if len(pdb_id) == 4:
        print(f"Downloading {pdb_id}...")
        url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
        data = requests.get(url).content
        if "not found" in str(data):
            print(f"Error: pdbID {pdb_id} is not found")
            return None
        with open(f"{id_dir}/{pdb_id}.pdb", "wb") as f:
            f.write(data)
            f.close()
            return os.path.abspath(f.name)
    else:
        print(f"Downloading {pdb_id}...")
        url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v3.pdb"
        data = requests.get(url).content
        if "not exist" in str(data):
            print(f"Error: uniprotID {pdb_id} is not found")
            return None
        with open(f"{id_dir}/AF-{pdb_id}-F1-model_v3.pdb", "wb") as f:
            f.write(data)
            f.close()
            return os.path.abspath(f.name)


def remove_alreadyexist_workbench(workbench_dir: str) -> None:
    """Delete directory if a workbench with the same name exists

    Args:
        workbench_dir (str): Path of workbench directory
    """
    if os.path.isdir(workbench_dir):
        print(f"{workbench_dir} already exists")
        shutil.rmtree(workbench_dir)
        print(f"Removed {workbench_dir}.")


def convert_complex_to_monomer(id_dir: str, pdb_path: str, output_pdb_name: str) -> str:
    """Convert (downloaded) pdb file to monomer

    Args:
        id_dir (str): Dictionary of ID and directory path
        pdb_path (str): Path of the pdb file to be converted to monomer
        output_pdb_name (str): Name of th pdb file converted to monomer

    Returns:
        str: Parh to the pdb file converted to monomer
    """
    pymol2_session = pymol2.PyMOL()
    pymol2_session.start()
    pymol2_session.cmd.load(pdb_path)
    pymol2_session.cmd.save(
        os.path.join(id_dir, output_pdb_name), "chain A and not resn HOH", format="pdb"
    )
    pymol2_session.stop()
    return os.path.join(id_dir, output_pdb_name)


def get_uniprotid_from_pdbid(pdb_id: str, pretty: str = False) -> str:
    """Obtain the corresponding uniprot ID from the PDB ID

    Args:
        pdb_id (str): PDB ID
        pretty (str, optional): Whether to output the result in a pretty format. Defaults to False.

    Returns:
        str: uniprot ID
    """
    pdb_id = pdb_id.lower()
    full_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id}?pretty={str(pretty).lower()}"
    json_results = requests.get(full_url).json()
    uniprot_id = json_results[pdb_id]
    uniprot_id2 = uniprot_id["UniProt"]
    uniprot_id3 = list(uniprot_id2.keys())
    assert len(uniprot_id3) == 1
    uniprot_id4 = uniprot_id3[0]
    return uniprot_id4


def download_fasta(pdb_id: str, output_dir: str) -> str:
    """Download fasta file from specified PDB or Uniprot ID

    Args:
        pdb_id (str): PDB or Uniprot ID
        output_dir (str): Path of the directory to save the fasta file obtained

    Returns:
        str: Path of the fasta file downloaded
    """
    output_dir = path_to_abspath(output_dir)
    uniprot_id = get_uniprotid_from_pdbid(pdb_id)
    print(f"Downloading FASTA file of {pdb_id} = {uniprot_id}...")
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

    fasta_path = f"{output_dir}/{pdb_id}.fasta"
    fasta_path = path_to_abspath(fasta_path)

    data = requests.get(url).content

    if "not exist" in str(data):
        print(f"Error occurred while downloading fasta file: pdbID {pdb_id} is not found")
        pass
    else:
        with open(fasta_path, "wb") as f:
            f.write(data)
            f.close()
    return fasta_path


def change_align_code_in_fasta(
    fasta_path: str, align_codes: str, alignment_format: str="fasta"
) -> str:
    """Parse fasta file and change align code

    Args:
        fasta_path (str): Path of fasta file
        align_codes (str): Align code
        alignment_format (str): Alignment format of fasta file. Defaults to "fasta"

    Returns:
        str: Align code
    """
    fasta_path = path_to_abspath(fasta_path)
    for record in SeqIO.parse(fasta_path,
                              alignment_format):
        id = record.id
        desc = record.description
        seq = record.seq
    record.id = align_codes
    SeqIO.write(record, fasta_path, alignment_format)
    return align_codes


def modelling_missing_res(
    pdb_id: str, id_dir: str, output_dir: str, fasta_path: str, pdb_path: str
) -> str:
    """Constructing Missing Regions by modeller. Output the highest scoring results as pdb file

    Args:
        pdb_id (str): PDB or Uniplot ID
        id_dir (str): Path to ID directory
        output_dir (str): Path of the directory where files related to the work done by the modeller are stored
        fasta_path (str): Path of fasta file
        pdb_path (str): Path of pdb file

    Returns:
        str: Path of the output pdb file of modeller
    """

    output_dir = path_to_abspath(output_dir)
    fasta_path = path_to_abspath(fasta_path)
    pdb_path = path_to_abspath(pdb_path)

    code_fill = change_align_code_in_fasta(
        fasta_path=fasta_path, align_codes=f"{pdb_id}_fill")

    modeller.log.none()
    env = modeller.Environ()
    aln = modeller.Alignment(env)
    env.io.atom_files_directory = [output_dir]
    aln.append(file=fasta_path, align_codes=code_fill, alignment_format="fasta")
    mdl = modeller.Model(env, file=pdb_path)
    aln.append_model(mdl, align_codes=pdb_id)
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False)
    aln.write(file=f"{output_dir}/{pdb_id}.ali", alignment_format="pir")

    os.chdir(output_dir)
    env = modeller.Environ()
    env.io.atom_files_directory = [id_dir]
    env.io.hetatm = True
    env.io.water = True
    aln_file = f"{output_dir}/{pdb_id}.ali"
    a = automodel.AutoModel(env, alnfile=aln_file, knowns=pdb_id, sequence=code_fill)
    a.md_level = automodel.refine.fast
    a.starting_model = 1
    a.ending_model = 2
    a.make()

    molpdf_list = []
    for i in range(len(a.outputs)):
        molpdf_list.append(a.outputs[i]["molpdf"])
    min_molpdf = min(molpdf_list)
    min_molpdf_index = molpdf_list.index(min_molpdf)
    min_model = a.outputs[min_molpdf_index]["name"]

    pdb_path_modelled = f"{output_dir}/{min_model}"
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    return pdb_path_modelled


config_path = "./config.ini"
flags.DEFINE_string(
    name="c", default=config_path, help="config file path. default is ./config.ini"
)
FLAGS = flags.FLAGS
if __name__ == "__main__":
    FLAGS(sys.argv)
    config_path = FLAGS.c
config = configparser.ConfigParser(allow_no_value=True, strict=False, delimiters="=")
config.optionxform = str
config.read(config_path)
flg_rm_disorder: bool = config.getboolean("SETTINGS", "remove_disordered_residue")
flg_insert_res: bool = config.getboolean("SETTINGS", "insert_residue_from_templete")
flg_packinkg: bool = config.getboolean("SETTINGS", "rosetta_packing")
flg_insert_sub: bool = config.getboolean("SETTINGS", "insert_substrate_from_templete")
flg_new_workbench: bool = config.getboolean(
    "SETTINGS", "make_brandnew_workbench_if_existed"
)
flg_comp_to_mono: bool = config.getboolean("SETTINGS", "convert_complex_to_monomer")
flg_model_missing: bool = config.getboolean("SETTINGS", "modelling_missing_residue")

workbench_dir = os.path.join(
    config["PATH"]["distination_path"], config["SETTINGS"]["workbench_dir_name"]
)
workbench_dir = path_to_abspath(workbench_dir)

if flg_new_workbench == True:
    remove_alreadyexist_workbench(workbench_dir=workbench_dir)

id_list: list = [key.upper() for key in config["ID"]]
id_dirs: dict = make_id_dirs(
    id_list=id_list,
    dist_dir=config["PATH"]["distination_path"],
    dir_name=config["SETTINGS"]["workbench_dir_name"]
)
id_pdb_paths: dict = {}


for ID, dir in id_dirs.items():
    id_pdb_paths[ID] = download_pdb_files(pdb_id=ID, id_dir=dir)

not_found_pdb_ids = []
for ID, pdb_path in id_pdb_paths.items():
    if pdb_path is None:
        os.removedirs(name=id_dirs[ID])
        not_found_pdb_ids = []
        not_found_pdb_ids.append(ID)

if len(not_found_pdb_ids) != 0:
    for ID in not_found_pdb_ids:
        print(f"{ID} Eliminated because of not found pdb file.")
        id_list.remove(ID)
        id_dirs.pop(ID)
        id_pdb_paths.pop(ID)

if flg_comp_to_mono is True:
    for ID, pdb_path in id_pdb_paths.items():
        if len(ID) == 4:
            print(f"Converting complex {ID} to monomer...")
            id_pdb_paths[ID] = convert_complex_to_monomer(
                id_dir=id_dirs[ID],
                pdb_path=pdb_path,
                output_pdb_name=f"{ID}_monomer.pdb"
            )

if flg_model_missing is True:
    for ID, pdb_path in id_pdb_paths.items():
        if len(ID) == 4:
            modeller_dir = os.path.join(id_dirs[ID], "modeller")
            os.mkdir(modeller_dir)
            fasta_path = download_fasta(pdb_id=ID, output_dir=modeller_dir)
            print(f"Modelling missing residue of {ID}...")
            id_pdb_paths[ID] = modelling_missing_res(
                pdb_id=ID,
                id_dir=id_dirs[ID],
                output_dir=modeller_dir,
                fasta_path=fasta_path,
                pdb_path=pdb_path
            )

if flg_rm_disorder is True:
    for ID, pdb_path in id_pdb_paths.items():
        print(f"Removing disordered residues from {ID}...")
        remove_disordered_residues(
            pdb_path=pdb_path,
            output_pdb_name=f"{id_dirs[ID]}/{ID}_disordered_removed.pdb",
        )
        id_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_disordered_removed.pdb"
else:
    print("Skip removing disordered residues")

if flg_insert_res is True:
    for ID, pdb_path in id_pdb_paths.items():
        templete_residue_name = config["RESIDUES_NAME_IN_TEMPLETE"][
            "insert_residue_name"
        ]
        print(f"Inserting {templete_residue_name} into {ID} from templete...")
        insert_templete_residue(
            residue_name=templete_residue_name,
            output_pdb_name=f"{id_dirs[ID]}/{ID}_res_inserted.pdb",
            pdb_path=pdb_path,
            templete_pdb_path=config["PATH"]["templete_pdb_path"],
        )
        id_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_res_inserted.pdb"
else:
    print("Inserting residues from templete is skipped")

if flg_packinkg is True:
    for ID, pdb_path in id_pdb_paths.items():
        print(f"Packing {ID} for optimizing sidechains and hetero atoms...")
        id_pdb_paths[ID] = rosetta_packing_residues(
            pdb_path=pdb_path,
            output_pdb_dir=id_dirs[ID],
            output_pdb_name=f"{ID}_packed.pdb",
        )
else:
    print("Sidechain packing is skipped")

if flg_insert_sub is True:
    substrate_name = config["RESIDUES_NAME_IN_TEMPLETE"]["insert_substrate_name"]
    for ID, pdb_path in id_pdb_paths.items():
        print(f"Inserting {substrate_name} into {ID} from templete...")
        insert_templete_residue(
            residue_name=substrate_name,
            output_pdb_name=f"{id_dirs[ID]}/{ID}_res_sub_inserted.pdb",
            pdb_path=pdb_path,
            templete_pdb_path=config["PATH"]["templete_pdb_path"],
        )
        id_pdb_paths[ID] = f"{id_dirs[ID]}/{ID}_res_sub_inserted.pdb"
else:
    print("Inserting substrate from templete is skipped")

for ID, pdb_path in id_pdb_paths.items():
    print(f"Excuting preparemd.py for {ID}...")
    preparemd_cmd = preparemd_settings_to_list(
        config=config,
        pdb_path=pdb_path,
        distdir=id_dirs[ID],
        preparemd_path=config["PATH"]["preparemd_script_path"],
    )
    subprocess.run(preparemd_cmd)

shutil.copy("config.ini", workbench_dir)
print("Copied config.ini to workbench directory.")

print("Making a script file for submitting MD jobs in yayoi...")
make_qscript(workbench_dir=workbench_dir, id_dirs=id_dirs)

print("Making a script file for extract init pdb and trajectry file...")
make_initscript(workbench_dir=workbench_dir, id_dirs=id_dirs)

print("Coping mdanalyze.py to workbench directory...")
shutil.copy("mdanalyze.py", workbench_dir)
shutil.copy("md_multianalyze.py", workbench_dir)

print("Process terminated.")
