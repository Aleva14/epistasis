import argparse
import os
import subprocess
import shutil


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
            # raise
    else:
        print("Directory ", path, "created.")


def foldx(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir, mut_id):
    global AA

    #FoldX command for BuildModel option. individual_list.txt - file with list of mutations.
    ind_list = "individual_list" + mut_id + ".txt" 
    foldx_command = "foldx --command=BuildModel --pdb-dir=%s --pdb=%s --mutant-file=" + ind_list + " --output-dir=%s --numberOfRuns=1 --out-pdb=false"

    #Create symlink to rotabase.txt (necessary for FoldX)
    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(OSError)
            raise
    else:
        print("Symlink to rotabase.txt created")

    #Create file with mutation
    ind_list_txt = open(ind_list, "w")
    mutation = aa1 + chain + pos + aa2
    ind_list_txt.write(mutation + ";")
    ind_list_txt.close()

    #Fill in the command with current values
    command = foldx_command % (pdb_folder, pdb + ".pdb", cur_dir)
    print(command)

    #Run FoldX
    proc = subprocess.run(command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)

    #foldx_output - list of lines of FoldX stdout
    foldx_output = proc.stdout.split("\n")
    
    #Check if run was ok.
    if "Specified residue not found." in proc.stdout:
        print("Wrong residue")
        return '-', 'Wrong residue'
    elif "No pdbs for the run found at:" in proc.stdout:
        print("PDB not found")
        return '-', 'PDB'
    
    #If ok - return average ddG of all runs
    else:
        ddG_result = cur_dir + "/Average_" + pdb + ".fxout"

        print(ddG_result)
        with open(ddG_result, 'r') as f:
            for line in f:
                if line.startswith(pdb):
                    ddG = line.split()[2]
    return ddG, 'No'    


def imutant(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir):

    # Command for I-Mutant
    # imut_command = "python /usr/local/share/I-Mutant2.0.7/I-Mutant2.0.py -pdb Test/1cei.pdb Test/1cei.dssp _ 17 A"
    imut_command_ddg =  "python -O /usr/local/share/I-Mutant2.0.7/I-Mutant2.0.py -pdbv %s %s %s %s %s"
    # Fill in the command with current values
    pdb_path = os.path.join(pdb_folder, pdb + ".pdb")
    dssp_path = os.path.join(cur_dir, pdb + ".dssp")
    command_ddg = imut_command_ddg % (pdb_path, dssp_path, chain, pos, aa2)

    # Run DSSP
    dssp(pdb_path, dssp_path)
    # Run I-MUTANT to find ddG
    proc = subprocess.Popen(command_ddg,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
    try:
        imut_out, errs = proc.communicate(timeout=3600)
    except subprocess.TimeoutExpired:
        proc.kill()
        return '-'
    imut_out = imut_out.strip().split("\n")
    result = imut_out[11].split()
    ddg = str(-float(result[3]))
  
    return ddg
 

def dssp(pdb_path, dssp_path):
    #command = "dssp pdb_short/1arr.pdb pdb_short/1arr.dssp"
    command = "dssp %s %s"
    cur_command = command % (pdb_path, dssp_path)
    print(cur_command)
    proc = subprocess.run(cur_command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)

def maestro(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir):
    # command = "maestro /usr/local/share/MAESTRO/config.xml pdb_short/3ssi.pdb --bu --evalmut='V13.A{I}'"
    command = "maestro /usr/local/share/MAESTRO/config.xml %s --bu --evalmut='%s%s.%s{%s}'"
    pdb_path = os.path.join(pdb_folder, pdb + ".pdb")
    maestro_command = command % (pdb_path, aa1, pos, chain, aa2)
    print(maestro_command)
    proc = subprocess.Popen(maestro_command,
                            shell=True,
                            universal_newlines=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        maestro_out, errs = proc.communicate(timeout=3600)
    except subprocess.TimeoutExpired:
        proc.kill()
        return '-'
    if "ERROR" in errs:
        return '-'
    maestro_out = maestro_out.strip().split("\n")[2]
    maestro_out = maestro_out.split()
    maestro_ddG = maestro_out[5]
    maestro_conf = maestro_out[6]
    maestro_len = maestro_out[1]

    return maestro_ddG


def main():
    global AA
    #Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="PDB identificator, for example 1C90")
    parser.add_argument("pdb_folder", help="Folder with pdb structures")
    parser.add_argument("chain", help="Protein chain, for example A")
    parser.add_argument("pos", help="Number of residue")
    parser.add_argument("aa1", help="One-letter initial aminoacid, for example P")
    parser.add_argument("aa2", help="One-letter mutated aminoacid, for example S")
    parser.add_argument("-directory", help="Directory to save results")
    args = parser.parse_args()
    if args.directory:
        path = args.directory
        create_dir(path)
    else:
        path = os.getcwd()

    pdb = args.pdb.lower()
    pos = args.pos
    aa1 = args.aa1
    aa2 = args.aa2
    pdb_folder = os.path.abspath(args.pdb_folder)
    chain = args.chain
    print(pdb, path, chain, pos, aa1, aa2, pdb_folder)
    #Turn all paths into absolute paths
    path = os.path.abspath(path)

    #Change working directory into path 
    os.chdir(path)

    #Create symlink to rotabase.txt (necessary for FoldX)
    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(OSError)
            raise
    else:
        print("Symlink to rotabase.txt created")
    

    result_path = os.path.join(path, pdb + aa1 + pos + chain + aa2 + ".csv")
    res_csv = open(result_path, 'w')
    
    mut_id = pdb + aa1 + pos + chain + aa2
    # Create directory to save output files of FoldX
    cur_dir = os.path.join(path, mut_id)
    print(cur_dir)
    create_dir(cur_dir)
    os.chdir(cur_dir)

    res_fields = "PDB,CHAIN,RES,AA1,AA2,Error,FX_ddG,IM_ddG,M_ddG\n"
    res_template = "%s,%s,%s,%s,%s,%s,%s,%s,%s\n"
    res_csv.write(res_fields)
    
    fx_ddG = "-"
    imut_ddG = "-"
    maestro_ddG = "-"


    #FoldX
    fx_ddG, err = foldx(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir, mut_id)
    if err == "No":
    # I-MUTANT
        imut_ddG = imutant(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir)

    # MAESTRO
        maestro_ddG = maestro(pdb, chain, pos, aa1, aa2, pdb_folder, cur_dir)
            
    res = res_template % (pdb, chain, pos, aa1, aa2, err, fx_ddG, imut_ddG, maestro_ddG)
    res_csv.write(res)
    res_csv.close()
    shutil.rmtree(cur_dir)


if __name__ == "__main__":
    main()
