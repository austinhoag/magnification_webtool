# magnification_webtool

The purpose of this webtool is to allow a user to upload a list of coordinates at which they want to query the magnifications in my lens models. 

The webtool front end can be viewed here: http://hoag.physics.ucdavis.edu/magnification_webtool/

The HTML form for the front-end is here: templates/glass.html

The webtool is launched from the script: webtool_glass.py

When the user uploads a list of coordinates and submits them, the back-end calculates the magnification from my models which are stored on a server at UCLA.

The back-end calculations are performed in the scripts in src/

## Requirements

Tested on python 2.6 and 2.7 

flask -- can be installed with pip

werkzeug -- can be installed with pip

astropy -- can be installed with pip (or automatically downloaded with the astroconda anaconda environment) 

cosmolopy -- can be installed with pip

## Example usage

Upload the file: example_upload.txt to the webtool and submit to see the magnification and 68% confidence intervals