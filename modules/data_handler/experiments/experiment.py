from time import localtime
from ..base_instance import BaseInstance

class Experiment(BaseInstance):
    def setup(self):
        self.content["done"] = False
        self.content["interrupted"] = False
        super().setup()

    def action_id(self, action):
        return self.content["pipeline"][action]["id"]

    def log_action_status(self, action, status):
        self.content["pipeline"][action][status] = localtime()
        super().store()

    def add_download(self, action, path):
        self.content["pipeline"][action]["file"] = path
        super().store()

    def start_action(self, action, cached = False):
        self.log_action_status(action, "started")
        if cached:
            self.content["pipeline"][action]["cached"] = True
            self.complete_action(action)

    def complete_action(self, action):
        self.log_action_status(action, "completed")

    def mark_error(self, action, error):
        self.content["error"] = str(error)
        self.log_action_status(action, "error")

    def mark_interrupted(self, action):
        self.content["interrupted"] = True
        self.log_action_status(action, "interrupted")

    def mark_done(self, experiment_id):
        self.content["done"] = localtime()
        super().store()
