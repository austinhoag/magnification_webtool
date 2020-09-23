# magnification_webtool

The purpose of this webtool is to allow a user to upload a list of coordinates at which they want to query the magnifications in my lens models. 

The HTML form for the front-end is here: templates/glass.html

The webtool is launched from the script: webtool_glass.py

When the user uploads a list of coordinates and submits them, the back-end calculates the magnification from my models which are stored on a server at UCLA.

The back-end calculations are performed in the scripts in src/

## Requirements

Tested on python 2.6 and 2.7 

pip install flask, werkzeug, astropy, cosmolopy

## Example usage

Upload the file: example_upload.txt to the webtool and submit to see the magnification and 68% confidence intervals
