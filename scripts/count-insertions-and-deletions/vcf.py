# Returns array [chr, pos, type, len]
def process_vcf_line(line):
    return None

#     line_result = empty_line_result()
#     if not line.startswith("#"):
#         vcf_fields = fields = line.split("\t")
#         if len(vcf_fields) < 10:
#             return -1
#         chromosome = vcf_fields[0]
#         line_result["chromosome"] = chromosome
#         reference_bases = vcf_fields[3]
#         alternative_bases = vcf_fields[4]
#         # Process insertions
#         if len(alternative_bases) > len(reference_bases):
#             line_result["counts"]["insertions"] = 1
#             line_result["counts"]["inserted_bases"] = len(alternative_bases) - len(reference_bases)
#         # Process deletions
#         if len(reference_bases) > len(alternative_bases):
#             line_result["counts"]["deletions"] = 1
#             line_result["counts"]["deleted_bases"] = len(reference_bases) - len(alternative_bases)
#     return line_result
