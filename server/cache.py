import os, json, uuid, sys

DATA_DIRECTORY = "../data/"
INDEX_PATH = DATA_DIRECTORY + "index.json"

class Cache:
    def __init__(self):
        self.__setup()
        with open(INDEX_PATH) as index_file:
            self.index = json.load(index_file)

    def __setup(self):
        if not os.path.isdir(DATA_DIRECTORY):
            os.mkdir(DATA_DIRECTORY)
            with open(INDEX_PATH, "w") as index_file:
                index_file.write(json.dumps({}));

    def __cache(self, url):
        dataset_id = str(uuid.uuid4())
        os.mkdir(DATA_DIRECTORY + dataset_id)
        self.index[url] = dataset_id
        with open(INDEX_PATH, "w") as index_file:
            index_file.write(json.dumps(self.index))
        return dataset_id

    def get_experiment(self, params):
        # TODO look for alignment and other results once implemented
        experiment = { "dataset": None }
        if params["dataset"] in list(self.index.keys()):
            experiment["dataset"] = dataset_path(self.index[params["dataset"]])
        return experiment

    def create_dataset(self, url):
        dataset_id = self.__cache(url)
        return dataset_path(dataset_id)

def dataset_path(dataset_id):
    return DATA_DIRECTORY + dataset_id + "/" + "data.fastq"

def is_uuid(id):
    try:
        uuid.UUID(dataset_id)
        return True
    except ValueError:
        return False
