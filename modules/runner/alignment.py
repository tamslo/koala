import os
import time
import datetime
import modules.file_utils as file_utils
from services import services

aligners = {}
for service in services:
    if service.type == "aligner":
        aligners[service.id] = service

def align(docker_client, data_handler, experiment, action_names):
    aligner_id = experiment.get("pipeline")[action_names["ALIGNMENT"]]["id"]
    if not aligner_id in aligners:
        raise Exception("The aligner with ID '{}' is not registered".format(aligner_id))
    aligner = aligners[aligner_id]

    # Prepare alignment
    destination = data_handler.cache.create_path(
        experiment,
        action_names["ALIGNMENT"]
    )
    file_utils.create_directory(destination)
    runtime_log_path = os.path.join(destination, "Runtime.log")

    # Define genome index path and temp path (will be renamed if successful)
    genome_index_path = data_handler.genome_index_path(experiment, aligner_id)
    temp_genome_index_path = genome_index_path + ".running"

    # If neccessary, build genome index
    if not os.path.exists(genome_index_path):
        try:
            index_parameters = {
                "docker_client": docker_client,
                "destination": destination,
                "genome_index_path": temp_genome_index_path,
                "reference_id": experiment.get("reference"),
                "reference_path": data_handler.reference_path(experiment)
            }
            build_genome_index_start = time.time()
            aligner.build_genome_index(index_parameters)
            with open(runtime_log_path, "w") as runtime_log:
                runtime = str(datetime.timedelta(
                    seconds=time.time() - build_genome_index_start)
                )
                runtime_log.write("Index generation: {}\n".format(runtime))
        except:
            file_utils.delete(temp_genome_index_path)
            raise
        os.rename(temp_genome_index_path, genome_index_path)
    else:
        with open(runtime_log_path, "w") as runtime_log:
            runtime_log.write("Index already present\n")

    # Run alignment
    out_file_name = "Out.sam"
    alignment_parameters = {
        "docker_client": docker_client,
        "destination": destination,
        "genome_index_path": genome_index_path,
        "dataset": data_handler.datasets.select(experiment.get("dataset")),
        "out_file_name": out_file_name
    }
    run_start = time.time()
    aligner.align(alignment_parameters)
    with open(runtime_log_path, "a") as runtime_log:
        runtime = str(datetime.timedelta(seconds=time.time() - run_start))
        runtime_log.write("Alignment: {}\n".format(runtime))

    # Create BAM file from SAM file
    sam_path = destination + out_file_name
    bam_path = destination + "Out.bam"
    docker_client.run_and_write_logs(
        "gatk",
        "samtools view -Sb /{}".format(sam_path),
        bam_path
    )

    return destination
