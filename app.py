import os, json, time
from flask import Flask, render_template

app = Flask(__name__)

@app.route("/")
def render_main():
    return render_template("layout.html")

if __name__ == "__main__":
    app.env = "development"
    app.debug = True
    app.run()
