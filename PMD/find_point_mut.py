import sys


AA = ('Asn', 'Leu', 'Glu', 'Asp', 'Trp', 'Lys', 'Met', 'Cys', 'Pro', 'Gly', 'Phe',
      'Ala', 'Tyr', 'His', 'Arg', 'Thr', 'Ile', 'Val', 'Gln', 'Ser')


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


def is_point_mutation(line):
    line = line.split()
    # print(line)
    if len(line) == 3:
        aa1 = line[0] in AA
        aa2 = line[2].strip() in AA
        pos = line[1].isnumeric()
        if aa1 and aa2 and pos:
            return True
    return False


def main(argv):
    global AA
    database = open_database(argv)
    output = open_output(argv)

    line = database.readline()
    n = 0
    while line:
        if line.startswith('CHANGE '):
            line = line.split(None, 1)
            if '(' in line[1]:
                print(line)
                parenth_pos = line[1].find('(')
                line[1] = line[1][:parenth_pos].strip()
                print(line)
            if is_point_mutation(line[1]):
                line = "CHANGE-POINT   " + line[1]
                print("IS POINT", line)
            else:
                line = line[0] + '  ' + line[1]
        if line.startswith('CHANGE-POINT'):
            if '(designated' in line or '(homozygous' in line:
                print(line)
                parenth_pos = line.find('(')
                line = line[:parenth_pos].strip() + '\n'
                print(line)
        output.write(line)
        line = database.readline()

    print(n)
    database.close()
    output.close()


if __name__ == "__main__":
    main(sys.argv)
