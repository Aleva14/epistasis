import argparse
import os
import subprocess


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


def foldx_repair(pdb, path, pdb_folder):
    foldx_command = "foldx --command=RepairPDB --pdb-dir=%s --pdb=%s --output-dir=%s"
    command = foldx_command % (pdb_folder, pdb, path)
    print(command)
    proc = subprocess.run(command,
                          shell=True,
                          universal_newlines=True,
                          stdout=subprocess.PIPE)
    print(proc.stdout)
    foldx_output = proc.stdout.split("\n")
    ddG = []
    time = 0
    for line in foldx_output:
        if line.startswith("Total"):
            if "time" in line:
                line = line.split()
                time = float(line[-2])
            else:
                line = line.split()
                ddG.append(float(line[-1]))
    print(ddG)
    print(time)
    return time, ddG[-1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_folder", help="Folder with pdb structures")
    parser.add_argument("-directory", help="Directory to save results, \
                        otherwise directory has name PDB_Repair")
    parser.add_argument("-n", help="Number of Repair_PDB runs. Default is 10")
    args = parser.parse_args()
    if args.directory:
        path = args.directory
    else:
        path = "PDB_Repair"

    if args.n:
        n_rounds = args.n
    else:
        n_rounds = 10

    path = os.path.abspath(path)
    args.pdb_folder = os.path.abspath(args.pdb_folder)
    create_dir(path)
    print(path)
    os.chdir(path)
    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(OSError)
            raise
    else:
        print("Symlink to rotabase.txt created")

    pdbs = os.listdir(args.pdb_folder)
    pdbs = list(map(lambda x: x.lower(), pdbs))
    print(pdbs)
    print(n_rounds)
    energy_output = path + "/energies.csv"
    energy_output = open(energy_output, 'w')
    head = "PDB,"
    for i in range(1, n_rounds):
        head += ",time" + str(i) + "," + str(i)
    head += "\n"
    energy_output.write(head)
    print(head)
    for structure in pdbs:
        result = structure + ","
        for i in range(n_rounds):
            time, ddG = foldx_repair(structure, path, args.pdb_folder)
            result += ','.join([str(time), str(ddG)])
        result += "\n"
        print(result)
        energy_output.write(result)
    energy_output.close()
if __name__ == "__main__":
    main()
