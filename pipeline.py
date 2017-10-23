import argparse
import os


foldx_command = "foldx --command=PositionScan --pdb-dir=%s --pdb=%s --positions=%s --output-dir=%s"


def create_dir(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            print("Unable to create directory ", path)
            print(os.strerror(OSError))
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
    mut_param = {}
    mut_param["No"] = mut[0]
    mut_param["PDB"] = mut[1].lower() + ".pdb"
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
    print(mut_param)
    return mut_param


def foldx(mut_param, pdb_folder, foldx_path):
    global foldx_command
    output_dir = foldx_path
    print(output_dir)
    mutation = mut_param["AA1"] + mut_param["CHAIN"] + mut_param["POS"] + mut_param["AA2"]
    command = foldx_command % (pdb_folder, mut_param["PDB"], mutation, output_dir)
    print(command)
    os.system(command)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mutations", help="File with mutations, format:\
                                            No\tPDB\tMUTATION\tCHAIN")
    parser.add_argument("pdb_folder", help="Folder with pdb structures")
    parser.add_argument("-directory", help="directory to save results, \
                    otherwise directory has same name as file with mutations")
    args = parser.parse_args()
    print(args.mutations)
    if args.directory:
        path = args.directory
    else:
        path = "NewFolder"

    try:
        mutations = open(args.mutations, 'r')
    except IOError:
        print('Can not open ', args.mutations)
    else:
        print(args.mutations, 'has', len(mutations.readlines()), 'lines')
        mutations.seek(0)

    create_dir(path)
    foldx_path = path + "/FoldX"
    create_dir(foldx_path)
    eris_path = path + "/Eris"
    create_dir(eris_path)

    try:
        os.symlink("/usr/local/share/FoldX/rotabase.txt", "rotabase.txt")
    except OSError:
        if not os.path.exists("rotabase.txt"):
            print(os.strerror(OSError))
            raise
    else:
        print("Symlink to rotabase.txt created")

    mut = mutations.readline()
    mut_param = parse_mut(mut)
    result = []
    while mut:
        foldx(mut_param, args.pdb_folder, foldx_path)
        mut = mutations.readline()
        mut_param = parse_mut(mut)

    mutations.close()


if __name__ == "__main__":
    main()
