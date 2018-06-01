import os, json, uuid, shutil

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

    def __write_index(self):
        with open(self.index_path, "w") as index_file:
            index_file.write(json.dumps(self.index))

    def __cache(self, url):
        dataset_id = str(uuid.uuid4())
        os.mkdir(self.directory + dataset_id)
        self.index[url] = dataset_id
        self.__write_index()
        return dataset_id

    def lookup(self, experiment, action):
        if experiment["dataset"] in self.index:
            dataset_id = self.index[experiment["dataset"]]
            if action == "dataset":
                return self.dataset_path(dataset_id)
            else:
                directory_path = self.directory + dataset_id + "/" + experiment[action]
                if os.path.isdir(directory_path):
                    return directory_path + "/" + os.listdir(directory_path)[0]
        # Default value
        return False

    def create_dataset(self, url):
        dataset_id = self.__cache(url)
        return self.dataset_path(dataset_id)

    def dataset_path(self, dataset_id):
        return self.directory + dataset_id + "/" + "data.fastq"

    def create_path(self, experiment, action):
        dataset_id = self.index[experiment["dataset"]]
        path = self.directory + dataset_id + "/" + experiment[action]
        os.makedirs(path)
        return path

    def clean_up(self, action, experiment):
        if not experiment["dataset"] in self.index:
            return None

        dataset_folder = self.directory + self.index[experiment["dataset"]]
        file_path = dataset_folder + "/" + experiment[action]
        if action == "dataset" and os.path.isdir(dataset_folder):
            shutil.rmtree(dataset_folder)
            del self.index[experiment[action]]
            self.__write_index()
        elif action != "dataset" and os.path.exists(file_path):
            os.remove(self.directory + dataset_folder + "/" + file_id)

def is_uuid(id):
    try:
        uuid.UUID(dataset_id)
        return True
    except ValueError:
        return False
