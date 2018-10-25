from time import localtime
from .base_instance import BaseInstance

class Experiment(BaseInstance):
    def setup(self):
        self.content["status"] = self.constants["experiment"]["WAITING"]
        super().setup()

    def action_id(self, action):
        return self.content["pipeline"][action]["id"]

    def log_action_status(self, action, status):
        self.content["pipeline"][action][status] = localtime()
        super().store()

    def add_download(self, action, path):
        self.content["pipeline"][action]["directory"] = path
        super().store()

    def start_action(self, action, cached = False):
        self.log_action_status(action, "started")
        if cached:
            self.content["pipeline"][action]["cached"] = True
            self.complete_action(action)

    def complete_action(self, action):
        self.log_action_status(action, "completed")

    def mark_error(self, action, error):
        self.content["status"] = self.constants["experiment"]["ERROR"]
        self.content["pipeline"][action]["message"] = str(error)
        self.log_action_status(action, "error")

    def mark_done(self):
        self.content["status"] = self.constants["experiment"]["DONE"]
        super().store()

    def mark_running(self):
        self.content["status"] = self.constants["experiment"]["RUNNING"]
        super().store()

    def mark_waiting(self):
        self.content["status"] = self.constants["experiment"]["WAITING"]
        super().store()

    def get_input_directory(self, action_id):
        pipeline = self.content["pipeline"]
        step_ids = []
        step_directories = []
        for step in list(map(lambda tuple: tuple[1], list(pipeline.items()))):
            for tuple in list(step.items()):
                if tuple[0] == "id":
                    step_ids.append(tuple[1])
                if tuple[0] == "directory":
                    step_directories.append(tuple[1])
        current_index = step_ids.index(action_id)
        if current_index == 0:
            raise Exception("The first action does not have a preceding action")
        return step_directories[current_index - 1]
