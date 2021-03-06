from crypt import methods
import os
import time
from flask import Flask, render_template, jsonify,  request
import json
from flask_cors import CORS
from orf.ORF_ultra_short_version import findORFs

app = Flask(__name__)
CORS(app)

def format_server_time():
    server_time = time.localtime()
    return time.strftime("%I:%M:%S %p", server_time)

@app.route("/hello", methods=['POST'])
def hello():
    return "Hello World!"

@app.route('/orf_api', methods=['POST'])
def orf_api():
    d = json.loads(request.get_data(as_text=True))
    orfs = findORFs(d.get('dna'),int(d.get('min')),int(d.get('max')))
    response =  jsonify(orfs)
    return response

@app.route('/')
def index():
    context = { 'server_time': format_server_time() }
    return render_template('index.html', context=context)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))