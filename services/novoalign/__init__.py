import os
from ..base_aligner import BaseAligner
import modules.file_utils as file_utils

class NovoAlign(BaseAligner):
    def __fasta_input(self, parameters):
        dataset = parameters["dataset"]
        evaluation = dataset.get("evaluation")
        return evaluation and evaluation["type"] == "beers"

    def __annotation_path(self, parameters):
        annotation_base_path = parameters["annotation_base_path"]
        reference_id = parameters["reference_id"]
        annotation_path = annotation_base_path + reference_id + ".txt"
        return os.path.exists(annotation_path) and annotation_path

    def __make_transcriptome_parameters(self, parameters):
        return {
            # sequence lenght radius, set to read length - 4bp
            "rad": int(parameters["dataset"].get("readLength")) - 4,
            # reduces the maximal number of splices per gene to be faster
            "num": 60000,
            # maximal minutes to process each gene's splices before interrupting
            "min": 10 # default, for the sake of completeness
        }

    def __masked_genome_path(self, parameters):
        reference_base_path = parameters["reference_base_path"]
        reference_id = parameters["reference_id"]
        return reference_base_path + reference_id + "_masked"

    def __junctions_prefix(self, parameters):
        annotation_base_path = parameters["annotation_base_path"]
        reference_id = parameters["reference_id"]
        make_transcriptome_parameters = self.__make_transcriptome_parameters(parameters)
        return "{}{}Rad{}Num{}kMin{}".format(
            annotation_base_path,
            reference_id,
            make_transcriptome_parameters["rad"],
            int(make_transcriptome_parameters["num"] / 1000),
            make_transcriptome_parameters["min"]
        )

    def __known_junctions_path(self, parameters):
        return "{}Splices.fasta".format(self.__junctions_prefix(parameters))

    def __theoretical_junctions_path(self, parameters):
        return "{}Transcripts.fasta".format(self.__junctions_prefix(parameters))

    def prepare_indexing(self, parameters):
        destination = parameters["destination"]
        reference_base_path = parameters["reference_base_path"]
        reference_id = parameters["reference_id"]
        annotation_path = self.__annotation_path(parameters)

        if annotation_path:
            # Mask genome
            masked_genome_path = self.__masked_genome_path(parameters)
            if not os.path.exists(masked_genome_path):
                command = "java -jar /opt/useq/Apps/MaskExonsInFastaFiles -f {}{} -u {} -s {}".format(
                    reference_base_path,
                    reference_id,
                    annotation_path,
                    masked_genome_path
                )
                output_parameters = { "log_file_path": destination + "Masking.log" }
                self.run_docker(command, parameters, output_parameters)

            # Create annotations
            known_junctions_path = self.__known_junctions_path(parameters)
            theoretical_junctions_path = self.__theoretical_junctions_path(parameters)
            if not os.path.exists(known_junctions_path):
                make_transcriptome_parameters = self.__make_transcriptome_parameters(parameters)
                # -s skips subsequent occurrences of splices with the same coordinates, memory intensive
                command = "java -jar /opt/useq/Apps/MakeTranscriptome" \
                    " -f /{} -u /{} -r {} -n {} -m {} -s".format(
                        reference_base_path + reference_id,
                        annotation_path,
                        make_transcriptome_parameters["rad"],
                        make_transcriptome_parameters["num"],
                        make_transcriptome_parameters["min"]
                    )
                output_parameters = { "log_file_path": destination + "Annotation.log" }
                self.run_docker(command, parameters, output_parameters)
                file_utils.unzip(known_junctions_path + ".gz")
                file_utils.unzip(theoretical_junctions_path + ".gz")

    def build_index_command(self, parameters):
        genome_index_path = parameters["genome_index_path"]
        reference_id = parameters["reference_id"]
        reference_path = parameters["reference_path"]
        annotation_path = self.__annotation_path(parameters)
        command = "novoindex -n {} /{} ".format(reference_id, genome_index_path)
        if annotation_path:
            command += "/{} /{} /{}".format(
                self.__masked_genome_path(parameters),
                self.__known_junctions_path(parameters),
                self.__theoretical_junctions_path(parameters),
            )
        else:
            command += reference_path
        return command

    def alignment_command(self, parameters):
        dataset = parameters["dataset"]
        genome_index_path = parameters["genome_index_path"]
        command = "novoalign -o SAM -f"
        file_utils.validate_file_content(genome_index_path)
        for direction, specification in dataset.get("data").items():
            command += " /{}".format(specification["path"])
        command += " -d /{}".format(genome_index_path)
        command += " -r All 10" # report max. 10 alignments per read
        command += " -v 0 70 70 '[>]([^:]*)'" # group junction and exon sequences together
        if self.__fasta_input(parameters):
            command += " -F FA"
        return command

    def conclude_alignment(self, parameters, sam_file_path):
        file_utils.validate_file_content(sam_file_path)
        if self.__annotation_path(parameters):
            command = "./fix_coordinates.sh {}".format(sam_file_path)
            output_parameters = { "log_file_path": parameters["destination"] + "Coordinates.log" }
            self.run_docker(command, parameters, output_parameters)
        # if self.__fasta_input(parameters):
            # TODO: Add dummy qualities for SAM File if FASTA mode

    def genome_index_amendment(self, parameters):
        data_handler = parameters["data_handler"]
        experiment = parameters["experiment"]
        dataset = data_handler.datasets.select(experiment.get("dataset"))
        if self.__annotation_path(parameters):
            return "_" + dataset.get("readLength")
        else:
            return super().genome_index_amendment(parameters)

class NovoAlignIndelSensitive(NovoAlign):
    def alignment_command(self, parameters):
        has_license = os.path.exists("services/novoalign/assets/novoalign.lic")
        if not has_license:
            raise Exception("NovoAlign in indel mode can only be run with license")

        command = super().alignment_command(parameters)
        command += " -x 3"
        command += " --matchreward 3"
        command += " --softclip 50,30" # four color
        # command += " --softclip 35,0" # two color
        return command
