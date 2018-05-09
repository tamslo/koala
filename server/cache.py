import os, json, uuid, sys

class Cache:
    def __init__(self, data_directory):
        self.directory = data_directory + "cache/"
        self.index_path = self.directory + "index.json"
        self.__setup()
        with open(self.index_path) as index_file:
            self.index = json.load(index_file)

    def __setup(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory)
            with open(self.index_path, "w") as index_file:
                index_file.write(json.dumps({}));

    def __cache(self, url):
        dataset_id = str(uuid.uuid4())
        os.mkdir(self.directory + dataset_id)
        self.index[url] = dataset_id
        with open(self.index_path, "w") as index_file:
            index_file.write(json.dumps(self.index))
        return dataset_id

    def get_experiment_data(self, params):
        # TODO look for alignment and other results once implemented
        experiment = { "dataset": None }
        if params["dataset"] in list(self.index.keys()):
            experiment["dataset"] = self.dataset_path(self.index[params["dataset"]])
        return experiment

    def create_dataset(self, url):
        dataset_id = self.__cache(url)
        return self.dataset_path(dataset_id)

    def dataset_path(self, dataset_id):
        return self.directory + dataset_id + "/" + "data.fastq"

def is_uuid(id):
    try:
        uuid.UUID(dataset_id)
        return True
    except ValueError:
        return False