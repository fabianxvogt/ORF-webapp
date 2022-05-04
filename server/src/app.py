from crypt import methods
import os
import time
from flask import Flask, render_template, jsonify,  request
from flask_restful import Api, Resource
import json
from flask_cors import CORS
from orf.ORF_ultra_short_version import findORFs

app = Flask(__name__)
CORS(app)
#api = Api(app)

def format_server_time():
    server_time = time.localtime()
    return time.strftime("%I:%M:%S %p", server_time)

@app.route("/hello", methods=['POST'])
def hello():
    return "Hello, Welcome to GeeksForGeeks"

@app.route('/orf_api', methods=['POST'])
def orf_api():
    #d = request.form

    d = json.loads(request.get_data(as_text=True))
    orfs = findORFs(d.get('dna'),int(d.get('min')),int(d.get('max')))
    response =  jsonify(orfs)
    #response = request.get_data(as_text=True)
    #response = jsonify({'some': 'data'})
    #response.headers.add('Access-Control-Allow-Origin', '*')
    return response

@app.route('/')
def index():
    context = { 'server_time': format_server_time() }
    return render_template('index.html', context=context)

if __name__ == '__main__':
    #api.add_resource(ORF_API, '/orf_api')  # '/users' is our entry point
    app.run(debug=True, host='0.0.0.0', port=int(os.environ.get('PORT', 8080)))