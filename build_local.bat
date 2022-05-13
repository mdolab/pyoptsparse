set CC=cl
set FC=flang
meson setup build --prefix=%cd%\installdir && meson install -C build