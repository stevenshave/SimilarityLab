import os
from config import Config
from flask import Flask, render_template, url_for, send_from_directory, redirect
from forms import *

app=Flask(__name__)
app.config.from_object(Config)
print(app.config['SECRET_KEY'])
@app.route('/favicon.ico')
def favicon():
    return send_from_directory(os.path.join(app.root_path, 'static'),
                               'favicon.ico', mimetype='image/vnd.microsoft.icon')

@app.route('/')
@app.route('/index.html')
def index():
    return render_template("index.html")

@app.route('/find_similars.html',  methods=['GET', 'POST'])
def find_similars():
    form=FindSimilarsForm()
    if form.validate_on_submit():
        return pos("/index.html",)
    return render_template("find_similars.html", form=form)

@app.route('/predict_targets.html')
def predict_targets():
    return render_template("predict_targets.html")

@app.route('/explore_compound.html')
def explore_compound():
    return render_template("explore_compound.html")

if __name__ == "__main__":
    app.run(debug=True)
