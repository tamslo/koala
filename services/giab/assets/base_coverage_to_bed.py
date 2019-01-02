import sys
import os

coverage_path = sys.argv[1]
min_coverage = int(sys.argv[2])
bed_path = sys.argv[3]
intermediate_path = coverage_path + ".{}-filtered".format(str(min_coverage))

if os.path.exists(intermediate_path):
    print("Filtered coverage file already present")
else:
    coverage_file = open(coverage_path, "r")
    intermediate_file = open(intermediate_path, "w")
    print("Writing filtered coverage file")
    for line in coverage_file:
        fields = line.split("\t")
        chromosome = fields[0]
        position = fields[1]
        try:
            coverage = int(fields[2])
        except Exception as exception:
            print("Error: {}".format(str(exception)))
            coverage = min_coverage
        finally:
            if coverage >= min_coverage:
                intermediate_file.write("{}\t{}\n".format(chromosome, position))
    coverage_file.close()
    intermediate_file.close()

intermediate_file = open(intermediate_path, "r")
bed_file = open(bed_path, "w")
print("Writing filtered BED file")

def get_fields(line):
    fields = line.split("\t")
    chromosome = fields[0]
    position = int(fields[1])
    return chromosome, position

def write_bed(file, chromosome, start_position, end_position):
    file.write("{}\t{}\t{}\n".format(chromosome, str(start_position), str(end_position)))

chromosome = None
start_position = None
end_position = None
first_line = True

for line in intermediate_file:
    if first_line:
        chromosome, start_position = get_fields(line)
        end_position = start_position
        first_line = False
        continue

    next_chromosome, next_position = get_fields(line)
    if next_chromosome == chromosome and (next_position - 1) == end_position:
        end_position = next_position
    else:
        write_bed(bed_file, chromosome, start_position, end_position)
        chromosome = next_chromosome
        start_position = next_position
        end_position = next_position

write_bed(bed_file, chromosome, start_position, next_position)

intermediate_file.close()
bed_file.close()
