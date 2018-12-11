import utils

comment_character = "#"
min_fields = 10

# Returns array of array [chr, pos, type, len]
def process_vcf_line(line):
    if not line.startswith(comment_character):
        fields = utils.get_fields(line, min_fields)
        if fields == None:
            return None
        chromosome = fields[0]
        position = int(fields[1])
        type = None
        length = None

        reference_bases = fields[3]
        alternative_bases = fields[4]
        if len(alternative_bases) > len(reference_bases): # Insertion
            type = utils.insertion_type
            length = len(alternative_bases) - len(reference_bases)
        elif len(reference_bases) > len(alternative_bases): # Deletion
            type = utils.deletion_type
            length = len(reference_bases) - len(alternative_bases)
        else:
            return None
        return [[chromosome, position, type, length]] # for compatibility
