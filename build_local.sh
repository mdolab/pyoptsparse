#!/usr/bin/bash
FC=gfortran meson setup build --prefix="$PWD/installdir"
meson install -C build