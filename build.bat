set CC=cl
set FC=flang
call python -c "import os, shutil; installpath = str(os.path.split(shutil.which('python'))[0]); os.system(f'meson setup build --prefix={installpath} && meson install -C build')"