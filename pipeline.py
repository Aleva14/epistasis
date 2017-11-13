import argparse
import os
import subprocess
import prody
import fnmatch

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
        return None, None
    mut = mut.split(',')
    if len(mut) < 5:
        print("Line has wrong format")
        return None, None
    mut_param = {}
    mut_param["No"] = mut[0]
    mut_param["PDB"] = mut[1].lower()
    mutation = mut[2].split(' ')
    if len(mutation) != 3:
        #print("Wrong mutation format for No ", mut_param["No"])
        return None, mut_param["No"]
    mut_param["AA1"] = mutation[0]
    mut_param["POS"] = mutation[1]
    mut_param["AA2"] = mutation[2]
    if (mut_param["AA1"] not in AA.keys()) or (mut_param["AA2"] not in AA.keys()):
        #print("Wrong mutation format for No ", mut_param["No"])
        return None, mut_param["No"]
    
    if mut[3] == '_':
        mut_param["CHAIN"] = 'A'
    else:
        mut_param["CHAIN"] = mut[3]
    mut_param["ddG_ProTherm"] = mut[4]
    # print(mut_param) 
    if (mut_param["AA1"] not in AA.keys()) or (mut_param["AA2"] not in AA.keys()):
        return None, "Wrong format"
    return mut_param, mut_param["No"]


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
        eris_out, errs = proc.communicate(timeout=3600)
    except subprocess.TimeoutExpired:
        proc.kill()
        return "timeout", '-', '-'
    if "Segmentation" in errs:
        return "segfault", '-', '-'
    if "design table" in eris_out:
        return "design", '-', '-' 

    #eris_out - list of lines of Eris stdout
    eris_out = eris_out.strip().split("\n")
    time = ''
    errs = errs.strip().split("\n")
    for line in errs:
        if line.startswith("real"):
            time = line.split()[-1]
    print(eris_out[-5].split()[-2])
    
    return "Ok", time, eris_out[-1].split()[-1]

def imutant(mut_param, pdb_folder):

    #Command for I-Mutant
    #imut_command = "python /usr/local/share/I-Mutant2.0.7/I-Mutant2.0.py -pdb Test/1cei.pdb Test/1cei.dssp _ 17 A"
    imut_command_sign = "python -O /usr/local/share/I-Mutant2.0.7/I-Mutant2.0.py -pdb %s %s %s %s %s"    
    imut_command_ddg =  "python -O /usr/local/share/I-Mutant2.0.7/I-Mutant2.0.py -pdbv %s %s %s %s %s"
    #Fill in the command with current values
    pdb_path = os.path.join(pdb_folder, mut_param["PDB"] + ".pdb")
    dssp_path = os.path.join(pdb_folder, mut_param["PDB"] + ".dssp")
    command_sign = imut_command_sign % (pdb_path, dssp_path, mut_param["CHAIN"], mut_param["POS"], mut_param["AA2"])
    command_ddg = imut_command_ddg % (pdb_path, dssp_path, mut_param["CHAIN"], mut_param["POS"], mut_param["AA2"])

    print(command_sign)
    
    #Run I-MUTANT to find sign
    proc = subprocess.Popen("time " + command_sign,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    try:
        imut_out, errs = proc.communicate(timeout=3600)
    except subprocess.TimeoutExpired:
        proc.kill()
        return "timeout", '-', '-', '-'

    imut_out = imut_out.strip().split("\n")
    time = ''
    errs = errs.strip().split("\n")
    for line in errs:
        if line.startswith("real"):
            time = line.split()[-1]
    result = imut_out[11].split()
    if result[3] == "Increase":
        sign = "-"
    else:
        sign = "+"
    RI = result[4]
    print("SIGN", sign)
    
    #Run I-MUTANT to find ddG
    proc = subprocess.Popen("time " + command_ddg,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    try:
        imut_out, errs = proc.communicate(timeout=3600)
    except subprocess.TimeoutExpired:
        proc.kill()
        return "timeout", '-', '-', '-'
    imut_out = imut_out.strip().split("\n")
    time = ''
    errs = errs.strip().split("\n")
    for line in errs:
        if line.startswith("real"):
            time = line.split()[-1]
    result = imut_out[11].split()
    ddg = str(-float(result[3]))
  
    return "Ok", time, sign, ddg, RI
 

def dssp(pdb_folder):
    #command = "dssp pdb_short/1arr.pdb pdb_short/1arr.dssp"
    command = "dssp %s %s"
    file_names = os.listdir(pdb_folder)
    for file in file_names:
        if fnmatch.fnmatch(file, "*.pdb"):
            dssp_name = file[:-3] + "dssp"
            if dssp_name not in file_names:
                cur_command = command % (os.path.join(pdb_folder, file), os.path.join(pdb_folder, dssp_name))
                proc = subprocess.run(cur_command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)


def main():
    global AA
    #Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("mutations", help="File with mutations, format:\
                                            No\tPDB\tMUTATION\tCHAIN")
    parser.add_argument("pdb_folder", help="Folder with pdb structures")
    parser.add_argument("-directory", help="directory to save results, \
                    otherwise directory has same name as file with mutations")
    parser.add_argument("-mut_nums", help="Numbers of mutations in the format \"1-5,8\"")
    args = parser.parse_args()
    if args.directory:
        path = args.directory
    else:
        path = "Experiment"

    if args.mut_nums:
        mut_nums = []
        print(args.mut_nums)
        args.mut_nums = args.mut_nums.split(",")
        for m in args.mut_nums:
            print(m)
            if m.isdigit():
                mut_nums.append(int(m))
            else:
                m = m.split("-")
                mut_nums += range(int(m[0]), int(m[1] + 1))
        if not mut_nums:
            mut_nums = "All"
    else:
        mut_nums = "All"
    print("Mutation numbers:",  mut_nums)            
    
    #Turn all paths into absolute paths
    path = os.path.abspath(path)
    mutations_path = os.path.abspath(args.mutations)
    pdb_folder_path = os.path.abspath(args.pdb_folder) 

    #Create folder to save results
    create_dir(path)
    foldx_path = os.path.join(path, "FoldX")
    create_dir(foldx_path)
    eris_path = os.path.join(path, "Eris")
    create_dir(eris_path)

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

    #Make dssp files for I-Mutant
    dssp(pdb_folder_path)

    result_path = os.path.join(path, "result.tsv")
    main_loop(result_path, foldx_path, eris_path, pdb_folder_path, mutations_path, mut_nums)


def main_loop(result_path, foldx_path, eris_path, pdb_folder, mutations_path, mut_nums):
    #Open list of mutations
    try:
        mutations = open(mutations_path, 'r')
    except IOError:
        print('Can not open ', mutations_path)
        raise
    else:
        print("Opened ", mutations_path)

    #Create csv with results
    try:
        res_csv = open(result_path, 'w')
    except IOError:
        print("Can not create result.tsv")
        raise
    else:
        print("Created result.tsv")

    res_fields = mutations.readline().strip() + ",FX,FX_time,FX_ddG,Eris,Eris_time,Eris_ddG,iMut,iMut_time,iMut_sign,iMut_ddg,iMut_RI\n"
    res_template = "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
    res_csv.write(res_fields)

    #Main loop
    mut = mutations.readline()
    mut_param, mut_num = parse_mut(mut)
    
    while mut:
        fx_res, fx_time, fx_ddG = "-", "-", "-"
        er_res, er_time, er_ddG = "-", "-", "-"
        imut, imut_time, imut_RI, imut_ddG = "-", "-", "-", "-"
        if not mut_num:
            fx_res = "Not valid"

        elif mut_nums == "All" or (int(mut_num) in mut_nums):
            print("Mutation number:", mut_num)
        
        #FoldX
            fx_res, fx_time, fx_ddG = foldx(mut_param, pdb_folder, foldx_path)
            print(mut_param["No"], ": ", fx_res, fx_time, fx_ddG)
            if fx_res == "Ok":

        # Eris
                er_res, er_time, er_ddG = eris(mut_param, pdb_folder, eris_path)
                print(mut_param["No"], ": ", er_res, er_time, er_ddG)
        

        # I-MUTANT
                imut, imut_time, imut_sign, imut_ddG, imut_RI = imutant(mut_param, pdb_folder)
                print(imut, imut_time, imut_sign, imut_ddG, imut_RI)
                print("\n")      
            res = res_template % (fx_res, fx_time, fx_ddG, er_res, er_time, er_ddG, imut, imut_time, imut_sign, imut_ddG, imut_RI)
            res_csv.write(mut.strip() + "," + res)

        #Read next mutation
        mut = mutations.readline()
        mut_param, mut_num = parse_mut(mut)

    os.system("mv " + "eris_results.txt" + " " + eris_path)
    mutations.close()
    res_csv.close()


if __name__ == "__main__":
    main()
