#!/bin/bash
python -c "import os, shutil; installpath = os.path.join(os.path.split(shutil.which('python'))[0], '..'); os.system(f'meson setup build --prefix={installpath} && meson install -C build')"