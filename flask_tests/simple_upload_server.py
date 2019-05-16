from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/')
def index():
  return render_template('upload.html')

# @app.route('/', methods=['POST'])

# def my_form_post():
#     text = request.form['text']
#     processed_text = text.upper()
#     return processed_text

# @app.route('/my-link/')
# def my_link():
#   print 'I got clicked!'

#   return 'Click.'

# @app.route('/calc',)
# def calc_mag():
# 	print "Calculating Magnification"
# 	return 'Calculated'

if __name__ == '__main__':
  app.run(debug=True)