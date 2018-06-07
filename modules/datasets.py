import os, json, uuid, shutil
from modules.file_handler import FileHandler

class Datasets:
    def __init__(self, data_directory):
        self.directory = data_directory + "datasets/"
        self.file_handler = FileHandler(self.directory)
        self.index_path = self.directory + "index.json"
        self.__setup()
        with open(self.index_path) as index_file:
            self.index = json.load(index_file)

    def __setup(self):
        if not os.path.exists(self.index_path):
            with open(self.index_path, "w") as index_file:
                index_file.write(json.dumps({}));

    def __write_index(self):
        with open(self.index_path, "w") as index_file:
            index_file.write(json.dumps(self.index))

    def select(self, dataset_id):
        return self.index[dataset_id]

    def create(self, dataset):
        os.mkdir(self.directory + dataset["id"])
        self.index[dataset["id"]] = dataset
        self.__write_index()
        return dataset

    def delete(self, dataset_id):
        self.file_handler.delete(self.dataset_path(dataset_id))
        dataset = self.index.pop(dataset_id)
        return dataset

    def lookup(self, experiment, action):
        dataset_id = experiment["dataset"]
        data_directory = self.directory + dataset_id + "/"
        if os.path.isdir(data_directory):
            dataset_path = self.dataset_path(dataset_id)
            if action == "dataset" and os.path.exists(dataset_path):
                return dataset_path
            else:
                action_directory = data_directory + experiment[action]
                if os.path.isdir(action_directory):
                    return data_directory + "/" + os.listdir(data_directory)[0]
        # Default value
        return False

    def dataset_path(self, dataset_id):
        return self.directory + dataset_id + "/" + "data.fastq"

    def create_path(self, experiment, action):
        dataset_id = experiment["dataset"]
        data_directory = self.directory + dataset_id + "/"
        if action == "dataset":
            os.makedirs(data_directory, exist_ok=True)
            path = self.dataset_path(dataset_id)
        else:
            path = self.directory + dataset_id + "/" + experiment[action]
            os.makedirs(path)
        return path

    def get_datasets(self):
        return self.index

    def clean_up(self, action, experiment):
        dataset_id = experiment["dataset"]
        dataset_folder = self.directory + dataset_id + "/"
        action_folder = dataset_folder + experiment[action]
        delete_path = dataset_folder if action == "dataset" else action_folder
        if os.path.isdir(delete_path):
            shutil.rmtree(delete_path)

def is_uuid(id):
    try:
        uuid.UUID(dataset_id)
        return True
    except ValueError:
        return False
