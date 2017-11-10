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
    global AA 
    if not mut:
        return None
    mut = mut.split(',')
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
    if (mut_param["AA1"] not in AA.keys()) or (mut_param["AA2"] not in AA.keys()):
        print("Wrong mutation format for No ", mut_param["No"])
        return None
    
    if mut[3] == '_':
        mut_param["CHAIN"] = 'A'
    else:
        mut_param["CHAIN"] = mut[3]
    mut_param["ddG_ProTherm"] = mut[4]
    # print(mut_param) 
    if (mut_param["AA1"] not in AA.keys()) or (mut_param["AA2"] not in AA.keys()):
        return None
    return mut_param


def foldx(mut_param, pdb_folder, foldx_path):
    global AA

    #FoldX command for BuildModel option. individual_list.txt - file with list of mutations.
    foldx_command = "foldx --command=BuildModel --pdb-dir=%s --pdb=%s --mutant-file=individual_list.txt --output-dir=%s --numberOfRuns=1 --out-pdb=false"

    #Create directory to save output files of FoldX separately for each mutation
    cur_dir = os.path.join(foldx_path, mut_param["No"])
    create_dir(cur_dir)

    #Create file with mutation
    ind_list = open("individual_list.txt", "w")
    mutation = mut_param["AA1"] + mut_param["CHAIN"] + mut_param["POS"] + mut_param["AA2"]
    ind_list.write(mutation + ";")
    ind_list.close()

    #Fill in the command with current values
    print(cur_dir)
    command = foldx_command % (pdb_folder, mut_param["PDB"] + ".pdb", cur_dir)
    print(command)

    #Run FoldX
    proc = subprocess.run(command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)

    #foldx_output - list of lines of FoldX stdout
    foldx_output = proc.stdout.split("\n")
    #for line in foldx_output:
    #     print(line)
    
    #Check if run was ok.
    if "Specified residue not found." in proc.stdout:
        return "Wrong residue", '-', '-'
    elif "No pdbs for the run found at:" in proc.stdout:
        return "PDB not found", '-', '-'
    
    #If ok - return time and average ddG of all runs
    else:
        for line in foldx_output:
            if line.startswith("Total time"):
                line = line.split()
                time = line[-2]
        ddG_result = cur_dir + "/Average_" + mut_param["PDB"] + ".fxout"
        print(ddG_result)
        with open(ddG_result, 'r') as f:
            for line in f:
                if line.startswith(mut_param["PDB"]):
                    ddG = line.split()[2]
        return "Ok", time, ddG


def eris(mut_param, pdb_folder, eris_path):

    #Command for Eris, including chain relaxation
    eris_command = "bash run_eris.sh ddg -m %s %s -r -f"
    
    #Extract only one chain from pdb
    atoms = prody.parsePDB(pdb_folder + "/" + mut_param["PDB"] + ".pdb", chain=mut_param["CHAIN"])
    pdb_name = mut_param["PDB"] + mut_param["CHAIN"] + ".pdb"
    prody.writePDB(pdb_name, atoms)
    
    #Fill in the command with current values
    command = eris_command % (mut_param["AA1"] + mut_param["POS"] + mut_param["AA2"], pdb_name)
    print(command)
    print("Eris calculating...")

    #Run Eris
    proc = subprocess.Popen("time " + command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    try:
        eris_out, errs = proc.communicate(timeout=1)
    except subprocess.TimeoutExpired:
        proc.kill()
        return "timeout", '-', '-'
    #print(errs) 
    if "Segmentation" in errs:
        return "segfault", '-', '-'
    if "design table" in eris_out:
        return "design", '-', '-' 

    #eris_out - list of lines of Eris stdout
    eris_out = eris_out.strip().split("\n")
    time = ''
    n = 0
    errs = errs.strip().split("\n")
    for line in errs:
        print(n, line)
        n += 1
        if line.startswith("real"):
            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
            time = line.split()[-1]
    print(eris_out[-5].split()[-2])
    os.system("mv " + "eris_results.txt" + " " + eris_path)
    print('TIME  ', time)
    return "Ok", time, eris_out[-1].split()[-1]

def main():
    global AA
    #Parse command line arguments
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
    
    #Turn all paths into absolute paths
    path = os.path.abspath(path)
    args.mutations = os.path.abspath(args.mutations)
    args.pdb_folder = os.path.abspath(args.pdb_folder) 

    #Open list of mutations
    try:
        mutations = open(args.mutations, 'r')
    except IOError:
        print('Can not open ', args.mutations)
        raise
    else:
        print("Opened ", args.mutations)

    #Create folder to save results
    create_dir(path)
    foldx_path = os.path.join(path, "FoldX")
    create_dir(foldx_path)
    eris_path = os.path.join(path, "Eris")
    create_dir(eris_path)

    try:
        res_csv = open(os.path.join(path, "result.tsv"), 'w')
    except IOError:
        print("Can not create result.tsv")
        raise
    else:
        print("Created result.tsv")

    #Change working directory into /tmp 
    os.chdir("/tmp/mc-buglakova")

    #Create symlink to rotabase.txt (necessary for FoldX)
    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(OSError)
            raise
    else:
        print("Symlink to rotabase.txt created")
    
    #Create symlink to the script which runs Eris
    try:
        os.symlink("/usr/local/share/Eris/run_eris.sh", "run_eris.sh")
    except OSError:
        if not os.path.exists("run_eris.sh"):
            print(OSError)
            raise
    else:
        print("Symlink to run_eris.sh created")


    res_fields = mutations.readline().strip() + ",FX,FX_time,FX_ddG,Eris,Eris_time,Eris_ddG\n"
    res_template = "%s,%s,%s,%s,%s,%s\n"
    res_csv.write(res_fields)

    #Main loop
    mut = mutations.readline()
    mut_param = parse_mut(mut)
    
    while mut:
        fx_res, fx_time, fx_ddG = "-", "-", "-"
        er_res, er_time, er_ddG = "-", "-", "-"
 
        if not mut_param:
            fx_res = "Not valid"

        # FoldX
        else:
            fx_res, fx_time, fx_ddG = foldx(mut_param, args.pdb_folder, foldx_path)
            print(mut_param["No"], ": ", fx_res, fx_time, fx_ddG)

        # Eris
            er_res, er_time, er_ddG = eris(mut_param, args.pdb_folder, eris_path)
            print(mut_param["No"], ": ", er_res, er_time, er_ddG)
              
        res = res_template % (fx_res, fx_time, fx_ddG, er_res, er_time, er_ddG)
        res_csv.write(mut.strip() + "," + res)

        #Read next mutation
        mut = mutations.readline()
        mut_param = parse_mut(mut)

    mutations.close()
    res_csv.close()


if __name__ == "__main__":
    main()
