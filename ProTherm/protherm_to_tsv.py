import sys

properties = ("NO.", "PROTEIN",  "MUTATION", "MUTATED_CHAIN", "ddG", "ddG_H2O", "PDB_wild", "PDB_mutant",  "T", "pH",
        "SOURCE", "LENGTH", "MOL-WEIGHT", "PIR_ID", "SWISSPROT_ID", "E.C.NUMBER", "EC", "PMD.NO", "NO_MOLECULE",
        "SEC.STR.", "BUFFER_NAME", "BUFFER_CONC", "ION_NAME_1", "ION_CONC_1", "ION_NAME_2", "ION_CONC_2", "ION_NAME_3", "ION_CONC_3",
        "ADDITIVES", "PROTEIN_CONC", "MEASURE", "METHOD", "dG_H2O", "dG", "Tm", "dTm", "dHcal", "m", "Cm", "dCp", "STATE", 
        "REVERSIBILITY", "ACTIVITY", "ACTIVITY_Km", "ACTIVITY_Kcat", "ACTIVITY_Kd", "KEY_WORDS", "REFERENCE", "AUTHOR", "REMARKS",
        "RELATED_ENTRIES", "ASA", "dHvH")

def open_database(argv):
    print(argv)
    try:
        database = open(argv[1], 'r')
        return database
    except IndexError:
        print("Usage: python3 protherm_to_tsv.py input (output)")
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


def read_block(database):
    global properties
    cur_line = database.readline()
    while not cur_line.startswith("NO."):
        cur_line = database.readline()
        if not cur_line:
            return cur_line
    cur_line = cur_line.rstrip()
    cur_line = cur_line.split(None, 1)
    
    block = {property: "-" for property in properties}
    #print(block)
    block["NO."] = cur_line[1]
    print(cur_line[1])
    n = cur_line[1]

    cur_line = database.readline()
    previous_property = "NO."
    while cur_line and not cur_line.startswith("//"):
        #print(cur_line)
        if n == "5927":
            print(cur_line)
        if cur_line.startswith("*"):
            pass
        else:
            cur_line = cur_line.rstrip()
            cur_line = cur_line.split(None, 1)
            if len(cur_line) == 2:
                if cur_line[0] not in properties:
                    block[previous_property] = block[previous_property] + ' '.join(cur_line)
                else:
                    block[cur_line[0]] = cur_line[1]
                    previous_property = cur_line[0]

        cur_line = database.readline()
    return block


def parse_database(database, output):
    global properties
    s = "\t"
    output.write(s.join(properties) + "\n")
    block = read_block(database)
    while block:
        row = [block[p] for p in properties]
        output.write(s.join(row) + "\n")
        block = read_block(database)

def main(argv):
    database = open_database(argv)
    output = open_output(argv)
    
    parse_database(database, output)

    database.close()
    output.close()


if __name__ == "__main__":
    main(sys.argv)
