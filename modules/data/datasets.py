import os, json, uuid, shutil
import modules.file_utils as file_utils

class Datasets:
    def __init__(self, data_directory, constants):
        self.directory = data_directory + "datasets/"
        self.index_path = self.directory + "index.json"
        self.constants = constants.dataset
        self.__setup()
        with open(self.index_path) as index_file:
            self.index = json.load(index_file)

    def __setup(self):
        file_utils.create_directory(self.directory)
        if not os.path.exists(self.index_path):
            file_utils.write(json.dumps({}), self.index_path)

    def __write_index(self):
        with open(self.index_path, "w") as index_file:
            index_file.write(json.dumps(self.index))

    def select(self, dataset_id):
        return self.index[dataset_id]

    def create(self, dataset, files=None):
        directory = self.directory + dataset["id"] + "/"
        os.mkdir(directory)
        for file_key in dataset["content"]:
            file_path = directory + file_key + ".fastq"
            if files:
                file = files[file_key]
                name = file.filename
                file.save(file_path)
                dataset["content"][file_key]["name"] = name
            dataset["content"][file_key]["path"] = file_path

        self.index[dataset["id"]] = dataset
        self.__write_index()
        return dataset

    def lookup(self, experiment, action):
        dataset_id = experiment["pipeline"]["dataset"]["id"]
        data_directory = self.dataset_folder(dataset_id)
        if os.path.isdir(data_directory):
            if action == "dataset":
                dataset = self.index[dataset_id]
                if len(os.listdir(data_directory)) > 0:
                    return data_directory
            else:
                action_id = experiment["pipeline"][action]["id"]
                action_directory = data_directory + action_id
                if os.path.isdir(action_directory):
                    return data_directory + "/" + os.listdir(data_directory)[0]
        # Default value
        return False

    def dataset_folder(self, dataset_id):
        return self.directory + dataset_id + "/"

    def create_path(self, experiment, action):
        dataset_id = experiment["pipeline"]["dataset"]["id"]
        action_id = experiment["pipeline"][action]["id"]
        data_directory = self.dataset_folder(dataset_id)
        path = self.dataset_folder(dataset_id) + action_id
        os.makedirs(path)
        return path

    def get_datasets(self):
        return self.index

    def clean_up(self, action, experiment):
        return None
        # dataset_id = experiment["pipeline"]["dataset"]["id"]
        # dataset = self.select(dataset_id)
        # dataset_folder = self.dataset_folder(dataset_id)
        # action_folder = dataset_folder + experiment["pipeline"][action]["id"]
        # delete_path = dataset_folder if action == "dataset" else action_folder
        # if action == "dataset" and dataset["method"] != self.constants["URL"]:
        #     return None
        # if os.path.isdir(delete_path):
        #     shutil.rmtree(delete_path)
        #     if action == "dataset":
        #         self.index.pop(dataset_id, None)

def is_uuid(id):
    try:
        uuid.UUID(dataset_id)
        return True
    except ValueError:
        return False
