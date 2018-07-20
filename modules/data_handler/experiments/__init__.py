from ..base_handler import BaseHandler
from .experiment import Experiment

class Experiments(BaseHandler):
    def __init__(self, directory):
        super().__init__(directory)
        self.Instance = Experiment
