from time import localtime
from services import get_services
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

    def __get_preceding_step_id(self, action_id):
        pipeline = self.content["pipeline"]
        step_ids = []
        action_ids = []
        for step_id, step in pipeline.items():
            step_ids.append(step_id)
            for tuple in list(step.items()):
                if tuple[0] == "id":
                    action_ids.append(tuple[1])
        current_index = action_ids.index(action_id)
        if current_index == 0:
            raise Exception("The first action does not have a preceding action")
        return step_ids[current_index - 1]

    def get_preceding_service(self, action_id):
        services = get_services()
        preceding_step = self.__get_preceding_step_id(action_id)
        action_id = self.content["pipeline"][preceding_step]["id"]
        return next(service for service in services if service.id == action_id)

    def get_input_directory(self, action_id):
        preceding_step = self.__get_preceding_step_id(action_id)
        return self.content["pipeline"][preceding_step]["directory"]
