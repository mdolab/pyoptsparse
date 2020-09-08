# OptView - Dash!

## Please place the OptViewDash folder found in the zip file in the postprocessing folder in pyoptsparse (where the current OptView.py file is found). It holds:

```
    - an assets folder holding styling files
    - OptView_dash2.py, the dash python file
    - 3 example history files (make sure to use a hist file compatible with your python version)
        - as_hist.hst is py2 compatible
        - tp109.hst is py3 compatible
        - opt2.hst is py3 compatible
    - Requirements.txt with needed dependencies
```

## You must have the following installed:

```
	- dash, plotly
        - To install, run: pip install -r requirements.txt
	- The most recent pyoptsparse
        - check instructions in the pyopstarse docs for updating/installing this
```

## To run, use this command:

```
	python OptView_dash2.py "history file name"

    Multiple history file support: You may also run the following command to display multiple history file data. For each file, it's data will be displayed with an _A, _B, _C, etc., in the order that you presented the file names as arguments. 
        - Run the following command for mult hist file support:
        - python OptView_dash2.py "hist file name #1" "hist file name #2" "hist file name #3"
```

To view the dash app, you will have to manually open the server in your browser that is listed in the terminal after running the above command

## Features Instructions

```
    - Auto-Refresh: This follows the same functionality as the old optview, allowing you to see the changes of an optimization as it is running.
        - If you toggle this checklist button, it will cause the program to default update every 10 seconds, however you may modify this refresh rate using the input box underneath 
        - I suggest making sure to toggle off this button when you are done or the optimization is complete so it does not add lag
        *Note*: This feature also works with multiple history files/optimizations running!
```
