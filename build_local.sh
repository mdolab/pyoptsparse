#!/usr/bin/bash
CC=gcc FC=gfortran meson setup build --prefix="$PWD/installdir"
meson install -C build