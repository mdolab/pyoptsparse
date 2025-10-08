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
We use `ruff <https://github.com/astral-sh/ruff>`_ and `pre-commit <https://github.com/pre-commit/pre-commit/>`_ for linting and formatting.
Please find the instructions in `this page <https://mdolab-mach-aero.readthedocs-hosted.com/en/latest/machFramework/contribute.html#coding-style>`_.

Additionally, we use `isort <https://github.com/PyCQA/isort>`_ for import sorting.
Please install it by

.. prompt:: bash

    pip install isort

and fix the sorting of all files by running:

.. prompt:: bash

    isort .

.. warning::
    For a PR to be accepted, it must pass all GitHub checks, which include formatting and syntax checks with ``ruff`` and ``pre-commit`` and import sorting checks with ``isort``.

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
