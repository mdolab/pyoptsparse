#!/usr/bin/bash
meson setup build --prefix="$PWD/installdir"
meson install -C build