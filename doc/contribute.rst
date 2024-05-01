How to Contribute to pyOptSparse
================================
pyOptSparse is an open-source tool, thus we welcome users to submit additions or fixes to the code to make it better for everybody.

Issues
------
If you have an issue with pyOptSparse, a bug to report, or a feature to request, submit an issue on the GitHub repository.
This lets other users know about the issue.
If you are comfortable fixing the issue, please do so and submit a pull request.

Coding style
------------
We use `black <https://github.com/psf/black>`_ for formatting.
Please install it following its documentation, and run it with the command line argument ``-l 120``, for example at the project root directory:

.. prompt:: bash

    black . -l 120

This will automatically format all Python files.

We use `flake8 <https://flake8.pycqa.org/en/latest/>`_ for linting.
Please install it following its instructions, and run it at the project root with:

.. prompt:: bash

    flake8 .

The configuration file we use for flake8 is a combination of `this file <https://github.com/mdolab/.github/blob/master/.flake8>`__ and the one at the root of this repository.
If there are any PEP-8 violations, ``flake8`` will print out the nature of the violation.
We run continuous integration with these tools on all pull requests submitted.
For an easier workflow, we recommend integrating these tools with your code editor.

.. warning::
    For a PR to be accepted, it must pass all GitHub checks, which include both formatting checks with ``black`` and syntax checks with ``flake8``.

Documentation
-------------
When you add or modify code, make sure to provide relevant documentation that explains the new code.
This should be done in code via docstrings and comments, but also in the Sphinx documentation as well if you add a new feature or capability.
Look at the ``.rst`` files in the ``doc`` section of the repo.

To build documentation locally, go to the ``doc`` folder and type ``make html``.
Building the documentation requires ``sphinx`` and ``numpydoc``, as well as the Sphinx RTD theme.
To install these dependencies, type

.. prompt:: bash

    pip install sphinx numpydoc sphinx-rtd-theme sphinx_mdolab_theme

Testing
-------
When you add code or functionality, add tests that cover the new or modified code.
These may be units tests for individual components or regression tests for entire models that use the new functionality.
All the existing tests can be found under the ``test`` folder.

Pull requests
-------------
Finally, after adding or modifying code, and making sure the steps above are followed, submit a pull request via the GitHub interface.
This will automatically go through all of the tests in the repo to make sure everything is functioning properly.
The main developers of pyOptSparse will then merge in the request or provide feedback on how to improve the contribution.
