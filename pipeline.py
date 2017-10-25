import argparse
import os
import subprocess
import prody

AA = {"A": "ALA", "R": "AGR", "N": "ASN", "D": "ASP", "B": "ASX", "C": "CYS",
      "E": "GLU", "Q": "GLN", "Z": "GLX", "G": "GLY", "H": "HIS", "I": "ILE",
      "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
      "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def create_dir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            print("Unable to create directory ", path)
            print(OSError)
            raise
        else:
            print("The directory ", path, "already exists")
            raise
    else:
        print("Directory ", path, "created.")


def parse_mut(mut):
    if not mut:
        return None
    mut = mut.split('\t')
    if len(mut) < 5:
        print("Line has wrong format")
        return None
    mut_param = {}
    mut_param["No"] = mut[0]
    mut_param["PDB"] = mut[1].lower()
    mutation = mut[2].split(' ')
    if len(mutation) != 3:
        print("Wrong mutation format for No ", mut_param["No"])
        return None
    mut_param["AA1"] = mutation[0]
    mut_param["POS"] = mutation[1]
    mut_param["AA2"] = mutation[2]
    if mut[3] == '_':
        mut_param["CHAIN"] = 'A'
    else:
        mut_param["CHAIN"] = mut[3]
    mut_param["ddG_ProTherm"] = mut[4]
    # print(mut_param)
    return mut_param


def foldx(mut_param, pdb_folder, foldx_path):
    global AA
    foldx_energy_terms = ["BackHbond", "SideHbond", "Energy_VdW", "Electro",
                          "Energy_SolvP", "Energy_SolvH", "Energy_vdwclash",
                          "energy_torsion", "backbone_vdwclash", "Entropy_sidec",
                          "Entropy_mainc", "water bonds", "helix dipole",
                          "loop_entropy", "cis_bond", "disulfide",
                          "kn electrostatic", "partial covalent interactions",
                          "Energy_Ionisation", "Entropy Complex"]

    foldx_command = "foldx --command=PositionScan --pdb-dir=%s --pdb=%s --positions=%s"

    mut_pdb = "%s_%s.pdb" % (AA[mut_param["AA2"]] + mut_param["POS"], mut_param["PDB"])
    wild_pdb = "%s_%s.pdb" % (AA[mut_param["AA1"]] + mut_param["POS"], mut_param["PDB"])
    empty = "PS_%s.fxout" % mut_param["PDB"]
    ddG_result = "PS_%s_scanning_output.txt" % mut_param["PDB"]
    energy = "energies_%s_%s.txt" % (mut_param["POS"], mut_param["PDB"])
    # print(mut_pdb, wild_pdb, empty, energy)

    mutation = mut_param["AA1"] + mut_param["CHAIN"] + mut_param["POS"] + mut_param["AA2"]
    command = foldx_command % (pdb_folder, mut_param["PDB"] + ".pdb", mutation)
    print(command)
    proc = subprocess.run(command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)
    foldx_output = proc.stdout.split("\n")
    # for line in foldx_output:
    #     print(line)
    if "Specified residue not found." in foldx_output:
        os.system("rm " + empty)
        os.system("rm " + energy)
        return "Wrong residue", None, None
    elif "No pdbs for the run found at:" in foldx_output:
        return "PDB not found", None, None
    else:
        cur_dir = foldx_path + "/" + mut_param["No"]
        create_dir(cur_dir)
        os.system("rm " + empty)
        os.system("mv " + wild_pdb + " " + cur_dir)
        os.system("mv " + mut_pdb + " " + cur_dir)

        energy_output = open(cur_dir + "/energies.txt", 'w')
        energy_output.write("File\t" + "\t".join(foldx_energy_terms) + "\n")
        with open(energy, 'r') as f:
            for line in f.readlines():
                energy_output.write(line)
        os.system("rm " + energy)

        dG = []
        for line in foldx_output:
            if line.startswith("Total          = "):
                line = line.split()
                # print(line)
                dG.append(line[-1])
        with open(ddG_result, 'r') as f:
            ddG = f.readlines()[1].split()[1]
        os.system("rm " + ddG_result)
        return "Ok", dG, ddG


def eris(mut_param, pdb_folder):
    eris_command = "bash run_eris.sh ddg -m %s %s -r -f"
    atoms = prody.parsePDB(pdb_folder + "/" + mut_param["PDB"] + ".pdb", chain=mut_param["CHAIN"])
    pdb_name = mut_param["PDB"] + mut_param["CHAIN"] + ".pdb"
    prody.writePDB(pdb_name, atoms)
    command = eris_command % (mut_param["AA1"] + mut_param["POS"] + mut_param["AA2"], pdb_name)
    print(command)
    print("Eris calculating...")
    proc = subprocess.run(command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)
    eris_output = proc.stdout.strip().split("\n")
    print(eris_output[-1].split()[-1])
    os.system("rm " + pdb_name)
    return eris_output[-1].split()[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mutations", help="File with mutations, format:\
                                            No\tPDB\tMUTATION\tCHAIN")
    parser.add_argument("pdb_folder", help="Folder with pdb structures")
    parser.add_argument("-directory", help="directory to save results, \
                    otherwise directory has same name as file with mutations")
    args = parser.parse_args()
    if args.directory:
        path = args.directory
    else:
        path = "Experiment"

    try:
        mutations = open(args.mutations, 'r')
    except IOError:
        print('Can not open ', args.mutations)
        raise
    else:
        print("Opened ", args.mutations)

    create_dir(path)
    foldx_path = path + "/FoldX"
    create_dir(foldx_path)
    eris_path = path + "/Eris"
    create_dir(eris_path)

    try:
        res_tsv = open(path + "/result.tsv", 'w')
    except IOError:
        print("Can not create result.tsv")
        raise
    else:
        print("Created result.tsv")

    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(OSError)
            raise
    else:
        print("Symlink to rotabase.txt created")

    res_fields = mutations.readline().strip() + "\tFX" + "\tFX_dG_WT0" +\
                 "\tFX_dG_WT1" + "\tFX_dG_MUT" + "\tFX_ddG\n"
    res_tsv.write(res_fields)
    mut = mutations.readline()
    mut_param = parse_mut(mut)
    while mut:
        if not mut_param:
            res_tsv.write(mut.strip() + "\t" + "Not valid" + "\t-\t-\t-\t-\n")
            mut = mutations.readline()
            mut_param = parse_mut(mut)
            continue

        # FoldX
        foldx_res, dG, ddG = foldx(mut_param, args.pdb_folder, foldx_path)
        print(mut_param["No"], ": ", foldx_res)
        if (foldx_res == "Wrong residue") or (foldx_res == "PDB not found"):
            foldx_res += "\t-\t-\t-\t-\t-\n"
            res_tsv.write(mut.strip() + "\t" + foldx_res)
            mut = mutations.readline()
            mut_param = parse_mut(mut)
            continue

        foldx_res = foldx_res + "\t" + "\t".join(dG) + "\t" + ddG + "\n"

        # Eris
        eris(mut_param, args.pdb_folder)
        res_tsv.write(mut.strip() + "\t" + foldx_res)
        mut = mutations.readline()
        mut_param = parse_mut(mut)

    mutations.close()
    res_tsv.close()


if __name__ == "__main__":
    main()
