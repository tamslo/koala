import json, zipfile, time, os, optparse
from flask import Flask, request, send_file, send_from_directory
from flask_cors import CORS
from collections import OrderedDict
from modules.constants import Constants
from modules.data.datasets import Datasets
from modules.data.experiments import Experiments
from modules.runner import Runner
from modules.exporter import Exporter

app = Flask(__name__)
CORS(app)

with open("services.json", "r") as services_file:
    services = json.load(services_file, object_pairs_hook=OrderedDict)

data_directory = "data/"
root_directory = os.path.dirname(os.path.abspath(__file__))

constants = Constants(root_directory)
datasets = Datasets(data_directory, constants)
experiments = Experiments(data_directory)
runner = Runner(datasets, experiments, data_directory, constants)
exporter = Exporter(data_directory)

@app.route("/ping")
def ping():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

@app.route("/context", methods=["GET"])
def get_context():
    return json.dumps({
        "services": services,
        "experiments": experiments.all(),
        "datasets": datasets.get_datasets()
    })

@app.route("/dataset", methods=["POST"])
def dataset():
    params = json.loads(request.form["json"])
    files = request.files
    dataset = datasets.create(params, files)
    return json.dumps(dataset)

@app.route("/experiment", methods=["POST", "PUT", "DELETE"])
def experiment():
    if request.method == "POST":
        params = request.get_json()
        experiment = experiments.create(params)
        return json.dumps(experiment)
    elif request.method == "PUT":
        params = request.get_json()
        experiment = experiments.update(params)
        return json.dumps(experiment)
    else:
        experiment_id = request.args.get("id")
        return json.dumps(experiments.delete(experiment_id))

@app.route("/execute", methods=["GET"])
def data():
    action = request.args.get("action")
    experiment_id = request.args.get("experiment")
    experiment = runner.execute(action, experiment_id)
    return json.dumps(experiment)

@app.route("/done", methods=["GET"])
def done():
    experiment_id = request.args.get("experiment")
    return json.dumps(experiments.mark_done(experiment_id))

@app.route("/export", methods=["GET"])
def export():
    experiment_id = request.args.get("experiment")
    path = request.args.get("path")
    if path != None:
        export_file_path = path
        export_file_name = exporter.file_name(path)
    elif experiment_id != None:
        experiment = experiments.select(experiment_id)
        export_file_path, export_file_name = exporter.zip(experiment)
    return send_file(
        export_file_path,
        as_attachment=True,
        attachment_filename=export_file_name
    )


# Routes for client built with `npm run build`

@app.route("/")
def serve():
    return send_from_directory("client/build/", "index.html")

@app.route("/koala.ico")
def servefav():
    return send_from_directory("client/build/", "koala.ico")

@app.route("/index.css")
def servecss():
    return send_from_directory("client/build/", "index.css")

@app.route("/static/js/<path:path>")
def servejs(path):
    return send_from_directory("client/build/static/js/", path)

@app.route("/static/media/<path:path>")
def servemedia(path):
    return send_from_directory("client/build/static/media/", path)

def clean_up():
    exporter.clean_up()
    # In case of an interruption, clean up experiments and datasets
    # If an error occurred, the components already handled it
    for experiment_id, experiment in experiments.all().items():
        if not experiment["done"] and not experiment["error"]:
            for action, pipeline_step in experiment["pipeline"].items():
                print(pipeline_step)
                started = "started" in pipeline_step and pipeline_step["started"]
                completed = "completed" in pipeline_step and pipeline_step["completed"]
                if started and not completed:
                    experiments.mark_interrupted(experiment_id, action)
                    datasets.clean_up(action, experiment)

if __name__ == "__main__":
    try:
        parser = optparse.OptionParser()
        parser.add_option("-d", "--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)
        options, _ = parser.parse_args()
        app.run(host="0.0.0.0", port="5000", debug=bool(options.debug))
    finally:
        clean_up()
