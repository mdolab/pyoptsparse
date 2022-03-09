Post-processing
===============
There are three post-processing utilities that are provided with pyOptSparse.

- OptView is a GUI designed to quickly and interactively visualize optimization histories
- OptView-Dash is a `Dash <https://plotly.com/dash/>`_ implementation of OptView
- ``History`` is a Python class that can be used to read in the history file, and provide API for programmatically extracting data.

.. _optview:

OptView
-------
Requirements
~~~~~~~~~~~~
``OptView`` has the following dependency tree::

    matplotlib (OptView)
      \-backends
      | \-backend_tkagg
      |   \-FigureCanvasTkAgg (OptView)
      |   \-NavigationToolbar2TkAgg (OptView)
      \-pyplot (OptView)
    mpl_toolkits
      \-axes_grid1
      | \-host_subplot (OptView)
      \-axisartist (OptView)
    numpy (OptView)

For installation instructions, see :ref:`install_optview`.  
Although not necessary for most usage, the ``dill`` package is needed if you wish to save an editable version of the graph produced in ``OptView``.
``dill`` can be installed via ``pip`` in a terminal using

.. prompt:: bash

    pip install dill

``view_saved_figure.py`` can be used to reformat and view the saved figure.

Usage
~~~~~
``OptView`` can be run via terminal from any directory as

.. prompt:: bash

    optview histFile

Here, ``histFile`` is the name of the history file to be examined (default is ``opt_hist.hst``).


Additionally, you can open multiple history files in the same ``OptView`` instance by calling them via the command line

.. prompt:: bash

    optview histFile1 histFile2 histFile3

Each file's contents will be loaded into ``OptView`` with a flag appended to the end of each variable or function name corresponding to the history file.
The first one listed will have '_A' added to the name, the second will have '_B' added, etc.
There is currently no limit to the number of history files than can be loaded.

Optionally, if you want to save the generated figures, there is an optional argument:

.. prompt:: bash

    optview histFile --output ~/my_figures

``outputDirectory`` is the name of the desired output directory for saved images.
By default, the figure is saved to the directory where you invoked ``optview``.

Features
~~~~~~~~
``OptView`` has many options and features, including:

    * plotting multiple variables on a single plot
    * producing stacked plots
    * live searchable variable names
    * hovering plot labels
    * saving the figure to an image or pickling it for later formatting
    * refreshing the optimization history on the fly

Although some of these are self-explanatory, the layout and usage of ``OptView`` will be explained below.

GUI Layout
~~~~~~~~~~
The window is divided into two sections.
The top is the canvas where the figure and graphs will be produced, while the bottom grayed section contains user-selectable options.
Here, we will focus on the user options.

The selectable variables are contained on the left hand side of the options panel in scrollable listboxes.
You can select multiple items from the listboxes using the normal selection operators such as control and shift.
If a selected variable is an array, a third listbox should appear on the right hand side of the options panel,
allowing you to select specific sub-variables within the single array variable.

There are three main options when selecting how to produce the graph(s):

    * Shared axes - all selected variables are plotted on a single pair of axes
    * Multiple axes - each selected variable gets its own y-axis while all selected data shares an x-axis
    * Stacked plots - each variable gets its own individual plot and the set is stacked vertically

Most checkbox options should play well with any of these three main options,
though there are known issues with using the 'multiple axes' option and delta values or for displaying arrays.

There are seven checkbox options:

    * Absolute delta values - displays the absolute difference between one iteration's value and the previous
    * Log scale - sets the y-axis as a log scale
    * Min/max for arrays - only shows the minimum and maximum value of a variable for each iteration
    * Show all for arrays - plots all variables within an array
    * Show legend - reveals the legend for the plotted data
    * Show bounds - shows the variable bounds as dashed lines
    * Show 'major' iterations - a filter to remove the line search iterations from the plotting results; especially useful for SNOPT output

Additionally, four buttons allow control of the plot:

    * Refresh history - reloads the history file; used if checking on an optimization run on the fly
    * Save all figures - saves .png versions of a basic plot for each variable in the history file
    * Save figure - saves a .png and .pickle version of the current plot (the .pickle version can be reformatted afterwards)
    * Quit - exits the program

Lastly, there are some miscellaneous features:

    * A search box to cull the selectable variables
    * A font size slider to control the text size on the plot
    * Hoverable tooltips when the cursor is on a plot line
    * A variable called `actual_iteration_number` that gives a translation between history file iteration number and run file iteration number. This is especially useful for debugging specific steps of an optimization or comparing values across different histories.


OptView-Dash
------------
This is a Dash_ implementation of OptView, and has many of the same features offered by OptView.
For installation instructions, see :ref:`install_optview`.
To run, use this command:

.. prompt:: bash

    optview_dash <filename>

Similar to OptView, you can invoke it with multiple history files.
To view the dash app, you will have to manually open the server in your browser that is listed in the terminal after running the above command.

Auto-refresh: This follows the same functionality as OptView, allowing you to see the changes of an optimization as it is running.

-  If you toggle this checklist button, it will cause the program to default update every 10 seconds, however you may modify this refresh rate using the input box underneath
-  Make sure to toggle off this button when you are done or the optimization is complete so it does not add lag.
-  This feature also works with multiple history files/optimizations running!

Directly Accessing the History Object
-------------------------------------
The history file generated by pyOptSparse is just a SqliteDict object.
To extract the stored information in Python, first initialize a History object:

.. code-block:: python

    hist = History("path/to/opt_hist.hst", flag="r")

From here, various information can be extracted, using the various ``get_`` methods.
To extract iteration history, use the function ``getValues()``.
See the page :ref:`history` for a full description of the history file structure and the API.
