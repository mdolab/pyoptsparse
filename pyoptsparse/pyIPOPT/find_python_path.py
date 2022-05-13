import os, shutil


def _find_python_path():
    if 'PYTHON' not in os.environ.keys():
        print(os.path.split(os.path.abspath(shutil.which("python")))[0])
    else:
        print(os.path.split(os.environ['PYTHON'])[0])


if __name__ == "__main__":
    _find_python_path()