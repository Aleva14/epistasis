import sys

properties = ("ENTRY", "CHANGE-POINT",  "-POINT", "STRUCTURE", "STABILITY",
              "SOURCE", "TITLE", "CROSS-REFERENCE", "EXPRESSION-SYSTEM",
              "PROTEIN", "SEQUENCE", "AUTHORS", "MEDLINE", "JOURNAL", "TRANSPORT",
              "FUNCTION", "DISEASE", "-REPLACE", "COMMENT", "COMMENT-1", "COMMENT-2",
              "N-TERMINAL", "CHANGE-REPLACE", "-DELETE", "MUTANT-NAME"
              "CHANGE-FRAGMENT", "-STOP", "-EXTEND", "EXTEND", "INSERT",
              "-MODIFY", "CHANGE-INSERT", "CHANGE-MODIFY", "CHANGE-DELETE",
              "PEPTIDE-SEQUENCE", "PROTEIN-SEQUENCE", "STRUCTURE/STABILITY")
# necessary_keys = ("ENTRY", "CHANGE-POINT", "STABILITY", "FUNCTION", "STRUCTURE",
#                  "STRUCTURE/STABILITY", "STABILITY/STRUCTURE", "STABILITY/TRANSPORT",
#                  "TRANSPORT/STABILITY", "EXPRESSION/STABILITY",
#                  "AUTHORS", "JOURNAL", "MEDLINE", "TITLE", "CROSS-REFERENCE",
#                  "PROTEIN", "SOURCE", "N-TERMINAL", "PROTEIN-SEQUENCE",
#                  "EXPRESSION-SYSTEM")
#
necessary_keys = ("ENTRY", "CHANGE-POINT", "STABILITY", "FUNCTION", "STRUCTURE",
                  "AUTHORS", "JOURNAL", "MEDLINE", "TITLE", "CROSS-REFERENCE",
                  "PROTEIN", "SOURCE", "N-TERMINAL", "PROTEIN-SEQUENCE",
                  "EXPRESSION-SYSTEM")

stability_keys = ("STRUCTURE/STABILITY", "STABILITY/STRUCTURE",
                  "TRANSPORT/STABILITY", "STABILITY/TRANSPORT",
                  "EXPRESSION/STABILITY", "STABILITY/EXPRESSION",
                  "STABILITY")

half_keys = ("-POINT", "-DELETE", "-REPLACE", "-STOP")

delimiters = ('CROSS-REFERENCE', 'CHANGE-POINT', 'CHANGE-INSERT',
              'CHANGE-FRAME', 'CHANGE-DELETE', 'CHANGE-EXTEND',
              'CHANGE-MODIFY', 'CHANGE-FRAGMENT', 'CHANGE-STOP',
              'CHANGE-REPLACE')


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


def dump(experiment, output):
    global necessary_keys
    global stability_keys

    experiment = {key: value for (key, value) in experiment}

    # Get all information about stability into STABILITY field
    stability = ''
    for key in stability_keys:
        if key in experiment.keys():
            stability += ' ' + experiment[key]
    if stability:
        experiment['STABILITY'] = stability

    if 'CHANGE-POINT' in experiment.keys() and 'STABILITY' in experiment.keys():
        result = ''
        for key in necessary_keys:
            if key in experiment.keys():
                result += experiment[key] + ','
            else:
                result += '-' + ','
            print(result[:-1] + '\n')

        output.write(result[:-1] + '\n')


def split_into_experiments(block, output):
    global delimiters
    experiment = []  # stack
    head_read = 0
    for prop in block:
        print("head_read", head_read)
        if prop[0] in delimiters and head_read == 1:
            dump(experiment, output)
            print(experiment)
            while experiment and experiment[-1][0] not in delimiters:
                print("DELIMITER", prop[0])
                print(experiment)
                experiment.pop()
            experiment.pop()
            if prop[0] == 'CROSS-REFERENCE':
                head_read = 0
        elif prop[0] in delimiters and head_read == 0 and prop[0] != 'CROSS-REFERENCE':
            head_read = 1
        experiment.append(prop)
    dump(experiment, output)


def read_block(database, output):
    global properties
    global half_keys
    global delimiters
    cur_line = database.readline()
    while not cur_line.startswith("ENTRY"):
        cur_line = database.readline()
        if not cur_line:
            return None

    block = []
    while cur_line and not cur_line.startswith("///"):
        print(cur_line)
        cur_line = cur_line.replace(',', '')
        if cur_line.startswith(' '):
            half_key = False
            for key in half_keys:
                if key in cur_line:
                    half_key = True
                    print(half_key)
                    cur_line = cur_line.split(None, 1)
                    print(cur_line)
                    block[-1][1] += '; ' + cur_line[1].strip()
            if not half_key:
                if block[-1][0] in delimiters:
                    block[-1][1] += ';' + cur_line.strip()
                else:
                    block[-1][1] += ' ' + cur_line.strip()

        else:
            cur_line = cur_line.split(None, 1)
            if len(cur_line) < 2:
                cur_line.append('-')
            if cur_line[0] in properties:
                block.append([cur_line[0], cur_line[1].strip()])
            else:
                print(cur_line[0])
        cur_line = database.readline().rstrip()
    for i in block:
        print(i)
    split_into_experiments(block, output)
    return block


def parse_database(database, output):
    global properties
    global necessary_keys
    s = ","
    output.write(s.join(necessary_keys) + "\n")
    block = read_block(database, output)
    while block:
        block = read_block(database, output)


def main(argv):
    database = open_database(argv)
    output = open_output(argv)

    parse_database(database, output)

    database.close()
    output.close()


if __name__ == "__main__":
    main(sys.argv)
