import sys


def open_database(argv):
    print(argv)
    try:
        database = open(argv[1], 'r', errors='replace')
        return database
    except IndexError:
        print("Usage: python3 pmd_to_csv.py input (output)")
        sys.exit(1)
    except IOError:
        print("Could not read file: ", argv[1])
        sys.exit(1)


def open_output(argv):
    if len(argv) >= 3:
        output_name = argv[2]
    else:
        output_name = 'output'
    try:
        output = open(output_name, 'w')
        return output
    except IOError:
        print("Could not open file: ", output_name)
        sys.exit(1)


def main(argv):
    database = open_database(argv)
    output = open_output(argv)

    line = database.readline()
    n = 0
    while line:
        stability_experiment = 0
        point_mutation = 0
        if line.startswith('ENTRY'):
            block = []
            block.append(line)
            line = database.readline()
            while line and not line.startswith('///'):
                print(line)
                block.append(line)
                if 'STABILITY' in line:
                    stability_experiment += 1
                    print(stability_experiment)
                if 'CHANGE-POINT' in line:
                    point_mutation += 1
                print(stability_experiment)
                line = database.readline()
            block.append("///\n")
        if stability_experiment and point_mutation:
            n += 1
            for l in block:
                output.write(l)
            stability_experiment = 0
            point_mutation = 0
        line = database.readline()

    print(n)
    database.close()
    output.close()


if __name__ == "__main__":
    main(sys.argv)
