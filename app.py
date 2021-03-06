import json, zipfile, time, os, optparse, atexit, yaml
from flask import Flask, request, send_file, send_from_directory
from flask_cors import CORS
from collections import OrderedDict
from modules.data_handler import DataHandler
from modules.runner import Runner
from modules.exporter import Exporter
from scripts.evaluation import collect_evaluation_results
from services import get_services

app = Flask(__name__)
CORS(app)

# Create constants.json for same constants in back- and frontend
with open("constants.yml", "r") as constants_file:
    constants = yaml.load(constants_file)
root_directory = os.path.dirname(os.path.abspath(__file__))
frontend_path = os.path.join(root_directory, "client/src/constants.json")
with open(frontend_path, "w") as constants_file:
    constants_file.write(json.dumps(constants))

data_directory = "data/"
data_handler = DataHandler(data_directory)
runner = Runner(data_handler, data_directory, constants)
exporter = Exporter(data_directory)

def get_content(instances):
    content = OrderedDict()
    for id, instance in instances.items():
        content[id] = instance.content
    return content

@app.route("/ping")
def ping():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

@app.route("/context", methods=["GET"])
def get_context():
    return json.dumps({
        "references": data_handler.get_references(),
        "datasets": get_content(data_handler.datasets.all()),
        "services": list(map(lambda Service: Service.frontend_information(), get_services())),
        "experiments": get_content(data_handler.experiments.all())
    })

@app.route("/dataset", methods=["POST"])
def dataset():
    params = json.loads(request.form["json"])
    files = request.files
    params["files"] = files
    dataset = data_handler.datasets.create(params)
    return json.dumps(dataset.content)

@app.route("/experiment", methods=["GET", "POST", "PUT", "DELETE"])
def experiment():
    if request.method == "GET":
        experiment_id = request.args.get("id")
        experiment = data_handler.experiments.select(experiment_id)
        return json.dumps(experiment.content)
    elif request.method == "POST":
        params = request.get_json()
        unpacked_pipeline = {}
        for action_name, action in params["pipeline"].items():
            if action == "":
                continue
            if isinstance(action, list):
                for index, action_id in enumerate(action):
                    name = action_name + "_" + str(index)
                    unpacked_pipeline[name] = { "id": action_id }
            else:
                unpacked_pipeline[action_name] = { "id": action }
        params["pipeline"] = unpacked_pipeline
        experiments = []
        datasets = params.pop("datasets")
        experiment_name = params.pop("name")
        for dataset in datasets:
            params["dataset"] = dataset
            params["name"] = "{} ({})".format(
                experiment_name,
                data_handler.datasets.select(dataset).get("name")
            )
            params.pop("id", None)
            experiment = data_handler.experiments.create(params)
            experiments.append(experiment.content)
            runner.add_task(experiment.get("id"))
        return json.dumps(experiments)
    elif request.method == "PUT":
        content = request.get_json()
        experiment = data_handler.experiments.update(content)
        runner.add_task(experiment.get("id"))
        return json.dumps(experiment.content)
    else: # request.method == "DELETE"
        experiment_id = request.args.get("id")
        if runner.current_task == experiment_id:
            raise Exception("Task cannot be deleted, it is already running")
        experiment = data_handler.experiments.delete(experiment_id)
        runner.remove_task(experiment_id)
        return json.dumps(experiment.content)

@app.route("/running", methods=["GET"])
def running():
    experiment_id = runner.current_task or len(runner.tasks) > 0 and runner.tasks[0]
    if experiment_id:
        content = data_handler.experiments.select(experiment_id).content
    else:
        content = {}
    return json.dumps(content)

@app.route("/export", methods=["GET"])
def export():
    experiment_id = request.args.get("experiment")
    path = request.args.get("path")
    if path != None:
        export_file_path, export_file_name = exporter.zip_experiment_directory(path)
    elif experiment_id != None:
        experiment = data_handler.experiments.select(experiment_id)
        export_file_path, export_file_name = exporter.zip_experiment(experiment)
    return send_file(
        export_file_path,
        as_attachment=True,
        attachment_filename=export_file_name
    )

@app.route("/evaluation", methods=["GET"])
def evaluation():
    evaluation_path = collect_evaluation_results(data_directory, data_handler)
    if len(os.listdir(evaluation_path)) == 0:
        return "No evaluation results present"

    export_name = "Evaluation.zip"
    export_path = exporter.zip(evaluation_path, export_name)
    return send_file(
        export_path,
        as_attachment=True,
        attachment_filename=export_name
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


if __name__ == "__main__":
    try:
        parser = optparse.OptionParser()
        parser.add_option("-d", "--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)
        options, _ = parser.parse_args()
        app.run(host="0.0.0.0", port="5000", debug=bool(options.debug), threaded=True)
    finally:
        exporter.clean_up()
        data_handler.clean_up()
        runner.scheduler.shutdown()
