import re
import utils

comment_character = "@"
min_fields = 11

# Returns lengths of operators (until index) as an array
def operator_lengths(cigar_string, index = None):
    if index != None:
        cigar_string = cigar_string[:index]
    digits = list(map(
        lambda string: int(string),
        re.findall(r'\d+', cigar_string)
    ))
    return digits

def indel_position(cigar_string, index):
    position = 0
    lengths = operator_lengths(cigar_string, index - 1)
    for length in lengths:
        position += length
    return position

def indel_length(cigar_string, index):
    return operator_lengths(cigar_string, index)[-1]

def get_from_cigar(cigar_string, operator, chromosome, position):
    matches = []
    indices = [match.start() for match in re.finditer(operator, cigar_string)]
    for index in indices:
        matches.append([
            chromosome,
            position + indel_position(cigar_string, index),
            operator,
            indel_length(cigar_string, index)
        ])
    return matches

def get_insertions(cigar_string, chromosome, position):
    return get_from_cigar(cigar_string, utils.insertion_type, chromosome, position)

def get_deletions(cigar_string, chromosome, position):
    return get_from_cigar(cigar_string, utils.deletion_type, chromosome, position)

def get_next_missing_one(cigar_string):
    return re.search(r'([A-Z][A-Z])', cigar_string)

# Chars to upper case and ddd 1s if none are given
def parse_cigar(cigar_string):
    cigar_string = cigar_string.upper()
    if re.search(r'^([A-Z])', cigar_string) != None:
        cigar_string = "1{}".format(cigar_string)
    next_missing_one = get_next_missing_one(cigar_string)
    while next_missing_one != None:
        missing_position = next_missing_one.start() + 1
        cigar_string = "{}1{}".format(
            cigar_string[:missing_position],
            cigar_string[missing_position:]
        )
        next_missing_one = get_next_missing_one(cigar_string)
    return cigar_string

# Returns array of arrays [chr, pos, type, len]
def process_sam_line(line):
    if not line.startswith(comment_character):
        fields = utils.get_fields(line, min_fields)
        if fields == None:
            return None
        chromosome = fields[2]
        position = int(fields[3])
        cigar_string = parse_cigar(fields[5])
        insertions = get_insertions(cigar_string, chromosome, position)
        deletions = get_deletions(cigar_string, chromosome, position)
        return insertions + deletions
