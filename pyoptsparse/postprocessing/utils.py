# Standard Python modules
import os

# External modules
import pkg_resources

ASSET_PATH = pkg_resources.resource_filename("pyoptsparse", os.path.join("postprocessing", "assets"))
