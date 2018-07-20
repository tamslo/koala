import json, zipfile, time, os, optparse
from flask import Flask, request, send_file, send_from_directory
from flask_cors import CORS
from collections import OrderedDict
from modules.constants import Constants
from modules.data_handler import DataHandler
from modules.runner import Runner
from modules.exporter import Exporter

app = Flask(__name__)
CORS(app)

data_directory = "data/"
data_handler = DataHandler(data_directory)
constants = Constants(os.path.dirname(os.path.abspath(__file__)))
runner = Runner(data_handler, data_directory, constants)
exporter = Exporter(data_directory)

with open("services.json", "r") as services_file:
    services = json.load(services_file, object_pairs_hook=OrderedDict)

@app.route("/ping")
def ping():
    return time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())

@app.route("/context", methods=["GET"])
def get_context():
    def get_content(instances):
        content = {}
        for id, instance in instances.items():
            content[id] = instance.content
        return content

    return json.dumps({
        "services": services,
        "experiments": get_content(data_handler.experiments.all()),
        "datasets": get_content(data_handler.datasets.all())
    })

@app.route("/dataset", methods=["POST"])
def dataset():
    params = json.loads(request.form["json"])
    files = request.files
    dataset = data_handler.datasets.create(params, files)
    return json.dumps(dataset.content)

@app.route("/experiment", methods=["GET", "POST", "PUT", "DELETE"])
def experiment():
    if request.method == "GET":
        experiment_id = request.args.get("id")
        experiment = data_handler.experiments.select(experiment_id)
        return json.dumps(experiment.content)
    elif request.method == "POST":
        params = request.get_json()
        experiment = data_handler.experiments.create(params)
        return json.dumps(experiment.content)
    elif request.method == "PUT":
        content = request.get_json()
        experiment = data_handler.experiments.update(content)
        return json.dumps(experiment.content)
    else: # request.method == "DELETE"
        experiment_id = request.args.get("id")
        experiment = data_handler.experiments.delete(experiment_id)
        return json.dumps(experiment.content)

@app.route("/execute", methods=["GET"])
def data():
    experiment_id = request.args.get("experiment")
    experiment = runner.execute(experiment_id)
    return json.dumps(experiment.content)

@app.route("/export", methods=["GET"])
def export():
    experiment_id = request.args.get("experiment")
    path = request.args.get("path")
    if path != None:
        export_file_path = path
        export_file_name = exporter.file_name(path)
    elif experiment_id != None:
        experiment = data_handler.experiments.select(experiment_id)
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
    

if __name__ == "__main__":
    try:
        parser = optparse.OptionParser()
        parser.add_option("-d", "--debug", action="store_true", dest="debug", help=optparse.SUPPRESS_HELP)
        options, _ = parser.parse_args()
        app.run(host="0.0.0.0", port="5000", debug=bool(options.debug))
    finally:
        exporter.clean_up()
        data_handler.clean_up()
