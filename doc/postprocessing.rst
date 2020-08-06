Post-processing
===============
There are two post-processing utilities that are provided with pyOptSparse.

- OptView is a GUI designed to quickly and interactively visualize optimization histories
- ``History`` is a Python class that can be used to read in the history file, and provide API for programmatically extracting data.

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

If you are successfully running pyOptSparse, these packages are most likely
already installed.

Although not necessary for most usage, the ``dill`` package is needed
if you wish to save an editable version of the graph produced in ``OptView``.
``dill`` can be installed via ``pip`` in a terminal using::

    pip install dill

``view_saved_figure.py`` can be used to reformat and view the saved figure.

Usage
~~~~~
``OptView`` can be run via terminal from any directory as::

    optview histFile

Here, ``histFile`` is the name of the history file to be examined
(default is ``opt_hist.hst``).


Additionally, you can open multiple history files in the same ``OptView`` instance
by calling them via the command line::

    optview histFile1 histFile2 histFile3

Each file's contents will be loaded into ``OptView`` with a flag appended to the end
of each variable or function name corresponding to the history file. The first one
listed will have '_A' added to the name, the second will have '_B' added, etc.
There is currently no limit to the number of history files than can be loaded.

Optionally, if you want to save the generated figures, there is an optional argument::

    optview histFile --output ~/my_figures

``outputDirectory`` is the name of the desired output directory for saved images.
By default, the figure is saved to the directory where you invokved ``optview``.

Features
~~~~~~~~
``OptView`` has many options and features, including:

    * plotting multiple variables on a single plot
    * producing stacked plots
    * live searchable variable names
    * hovering plot labels
    * saving the figure to an image or pickling it for later formatting
    * refreshing the optimization history on the fly

Although some of these are self-explanatory, the layout and usage of ``OptView``
will be explained below.

GUI Layout
~~~~~~~~~~
The window is divided into two sections.
The top is the canvas where the figure and graphs will be produced,
while the bottom grayed section contains user-selectable options.
Here, we will focus on the user options.

The selectable variables are contained on the left hand
side of the options panel in scrollable listboxes.
You can select multiple items from the listboxes using the normal selection
operators such as control and shift.
If a selected variable is an array, a third listbox should appear on the
right hand side of the options panel, allowing you to select specific
sub-variables within the single array variable.

There are three main options when selecting how to produce the graph(s):

    * Shared axes - all selected variables are plotted on a single pair of axes
    * Multiple axes - each selected variable gets its own y-axis while all selected data shares an x-axis
    * Stacked plots - each variable gets its own individual plot and the set is stacked vertically

Most checkbox options should play well with any of these three main options,
though there are known issues with using the 'multiple axes'
option and delta values or for displaying arrays.

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

More features are being developed on an as-needed basis.
Feel free to edit the code as you see fit and submit a pull request if you
would like to see a feature added.
Alternatively, you can submit an issue ticket to discuss possible features.

Directly Accessing the History Object
-------------------------------------
The history file generated by pyOptSparse is just a SqliteDict object.
To extract the stored information in Python, first initialize a History object:

.. code-block:: python

    >>> hist = History('path/to/opt_hist.hst', flag='r')

From here, various information can be extracted, using the various ``get_`` methods.
To extract iteration history, use the function ``getValues()``.
See the page :ref:`history` for a full description of the history file structure and the API.
