from flask import Flask, render_template, request
from calc_magnification_glass import get_magnification_list_smart
from werkzeug import secure_filename

app = Flask(__name__)

@app.route('/')
def index():
  return render_template('glass.html')

@app.route('/', methods=['POST'])
def magnifications_from_file():
	try:
		input_file = request.files['upfile']
		filename = secure_filename(input_file.filename)
		savename = './uploads/%s' % filename
		input_file.save(savename)
		print "Saved %s" % filename
		return get_magnification_list_smart(savename)
	except:
		return "There was an error. Check coordinates"
if __name__ == '__main__':
  app.run(debug=True,threaded=True)
