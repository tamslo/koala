import os
import modules.file_utils as file_utils
import services.star.alignment as star
import services.novoalign.alignment as novoalign
import services.test_aligner.alignment as test_aligner

aligners = {
    "test_aligner": test_aligner,
    "star": star,
    "novoalign": novoalign
}

def align(docker_client, data_handler, experiment, action_names):
    destination = data_handler.cache.create_path(
        experiment,
        action_names["ALIGNMENT"]
    )
    aligner_id = experiment.get("pipeline")[action_names["ALIGNMENT"]]["id"]
    if not aligner_id in aligners:
        raise Exception("The aligner with ID '{}' is not registered".format(aligner_id))
    aligner = aligners[aligner_id]

    file_utils.create_directory(destination)
    genome_index_path = data_handler.genome_index_path(experiment, aligner_id)
    temp_genome_index_path = genome_index_path + ".running"

    if not os.path.exists(genome_index_path):
        try:
            aligner.build_genome_index(
                aligner_id,
                docker_client,
                temp_genome_index_path,
                data_handler,
                experiment,
                destination
            )
        except:
            file_utils.delete(temp_genome_index_path)
            raise
        os.rename(temp_genome_index_path, genome_index_path)

    dataset = data_handler.datasets.select(experiment.get("dataset"))
    aligner.run(aligner_id, docker_client, dataset, genome_index_path, destination)
    return destination
