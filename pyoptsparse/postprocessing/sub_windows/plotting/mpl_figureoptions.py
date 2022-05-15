# Copyright © 2009 Pierre Raybaut
# Licensed under the terms of the MIT License
# see the Matplotlib licenses directory for a copy of the license


"""Module that provides a GUI-based editor for Matplotlib's figure options."""

# External modules
# import matplotlib
from matplotlib import cbook, cm
from matplotlib import colors as mcolors
from matplotlib import image as mimage
from matplotlib import markers
from matplotlib.backends.qt_compat import QtGui
from matplotlib.backends.qt_editor import _formlayout
from matplotlib.dates import DateConverter, num2date

LINESTYLES = {
    "-": "Solid",
    "--": "Dashed",
    "-.": "DashDot",
    ":": "Dotted",
    "None": "None",
}

DRAWSTYLES = {
    "default": "Default",
    "steps-pre": "Steps (Pre)",
    "steps": "Steps (Pre)",
    "steps-mid": "Steps (Mid)",
    "steps-post": "Steps (Post)",
}

MARKERS = markers.MarkerStyle.markers


def figure_edit(axes, parent=None):
    """Edit matplotlib figure options"""
    sep = (None, None)  # separator

    # Get / General
    def convert_limits(lim, converter):
        """Convert axis limits for correct input editors."""
        if isinstance(converter, DateConverter):
            return map(num2date, lim)
        # Cast to builtin floats as they have nicer reprs.
        return map(float, lim)

    xconverter = axes.xaxis.converter
    xmin, xmax = convert_limits(axes.get_xlim(), xconverter)
    yconverter = axes.yaxis.converter
    ymin, ymax = convert_limits(axes.get_ylim(), yconverter)
    general = [
        ("Title", axes.get_title()),
        sep,
        (None, "<b>X-Axis</b>"),
        ("Left", xmin),
        ("Right", xmax),
        ("Label", axes.get_xlabel()),
        ("Scale", [axes.get_xscale(), "linear", "log", "symlog", "logit"]),
        sep,
        (None, "<b>Y-Axis</b>"),
        ("Bottom", ymin),
        ("Top", ymax),
        ("Label", axes.get_ylabel()),
        ("Scale", [axes.get_yscale(), "linear", "log", "symlog", "logit"]),
        sep,
        ("Show Legend", True if axes.legend_ is not None else False),
    ]

    if axes.legend_ is not None:
        old_legend = axes.get_legend()
        _draggable = old_legend._draggable is not None
        _ncol = old_legend._ncol
        _fontsize = int(old_legend._fontsize)
        _frameon = old_legend.get_frame_on()
        _framecolor = mcolors.to_hex(mcolors.to_rgba(old_legend.get_frame().get_edgecolor()), keep_alpha=True)
        _bgcolor = mcolors.to_hex(
            mcolors.to_rgba(old_legend.get_frame().get_facecolor(), old_legend.get_frame().get_alpha()),
            keep_alpha=True,
        )
        _framealpha = old_legend.get_frame().get_alpha()
    else:
        _draggable = False
        _ncol = 1
        _fontsize = 15
        _frameon = False
        _framecolor = "000000"
        _bgcolor = "#ffffff"
        _framealpha = 0.5

    legend = [
        ("Draggable", _draggable),
        ("columns", _ncol),
        ("Font Size", _fontsize),
        ("Frame", _frameon),
        ("Frame Color", _framecolor),
        ("Background Color", _bgcolor),
        ("Alpha", _framealpha),
    ]

    # Save the unit data
    xunits = axes.xaxis.get_units()
    yunits = axes.yaxis.get_units()

    # Get / Curves
    labeled_lines = []
    for line in axes.get_lines():
        label = line.get_label()
        if label == "_nolegend_":
            continue
        labeled_lines.append((label, line))
    curves = []

    def prepare_data(d, init):
        """
        Prepare entry for FormLayout.

        *d* is a mapping of shorthands to style names (a single style may
        have multiple shorthands, in particular the shorthands `None`,
        `"None"`, `"none"` and `""` are synonyms); *init* is one shorthand
        of the initial style.

        This function returns an list suitable for initializing a
        FormLayout combobox, namely `[initial_name, (shorthand,
        style_name), (shorthand, style_name), ...]`.
        """
        if init not in d:
            d = {**d, init: str(init)}
        # Drop duplicate shorthands from dict (by overwriting them during
        # the dict comprehension).
        name2short = {name: short for short, name in d.items()}
        # Convert back to {shorthand: name}.
        short2name = {short: name for name, short in name2short.items()}
        # Find the kept shorthand for the style specified by init.
        canonical_init = name2short[d[init]]
        # Sort by representation and prepend the initial value.
        return [canonical_init] + sorted(short2name.items(), key=lambda short_and_name: short_and_name[1])

    for label, line in labeled_lines:
        color = mcolors.to_hex(mcolors.to_rgba(line.get_color(), line.get_alpha()), keep_alpha=True)
        ec = mcolors.to_hex(mcolors.to_rgba(line.get_markeredgecolor(), line.get_alpha()), keep_alpha=True)
        fc = mcolors.to_hex(mcolors.to_rgba(line.get_markerfacecolor(), line.get_alpha()), keep_alpha=True)
        curvedata = [
            ("Label", label),
            sep,
            (None, "<b>Line</b>"),
            ("Line style", prepare_data(LINESTYLES, line.get_linestyle())),
            ("Draw style", prepare_data(DRAWSTYLES, line.get_drawstyle())),
            ("Width", line.get_linewidth()),
            ("Color (RGBA)", color),
            sep,
            (None, "<b>Marker</b>"),
            ("Style", prepare_data(MARKERS, line.get_marker())),
            ("Size", line.get_markersize()),
            ("Face color (RGBA)", fc),
            ("Edge color (RGBA)", ec),
        ]
        curves.append([curvedata, label, ""])
    # Is there a curve displayed?
    has_curve = bool(curves)

    # Get ScalarMappables.
    labeled_mappables = []
    for mappable in [*axes.images, *axes.collections]:
        label = mappable.get_label()
        if label == "_nolegend_" or mappable.get_array() is None:
            continue
        labeled_mappables.append((label, mappable))
    mappables = []
    cmaps = [(cmap, name) for name, cmap in sorted(cm._colormaps.items())]
    for label, mappable in labeled_mappables:
        cmap = mappable.get_cmap()
        if cmap not in cm._colormaps.values():
            cmaps = [(cmap, cmap.name), *cmaps]
        low, high = mappable.get_clim()
        mappabledata = [
            ("Label", label),
            ("Colormap", [cmap.name] + cmaps),
            ("Min. value", low),
            ("Max. value", high),
        ]
        if hasattr(mappable, "get_interpolation"):  # Images.
            interpolations = [(name, name) for name in sorted(mimage.interpolations_names)]
            mappabledata.append(("Interpolation", [mappable.get_interpolation(), *interpolations]))
        mappables.append([mappabledata, label, ""])
    # Is there a scalarmappable displayed?
    has_sm = bool(mappables)

    datalist = [(general, "Axes", ""), (legend, "Legend", "")]
    if curves:
        datalist.append((curves, "Curves", ""))
    if mappables:
        datalist.append((mappables, "Images, etc.", ""))

    def apply_callback(data):
        """A callback to apply changes."""
        orig_xlim = axes.get_xlim()
        orig_ylim = axes.get_ylim()

        general = data.pop(0)
        legend = data.pop(0)
        curves = data.pop(0) if has_curve else []
        mappables = data.pop(0) if has_sm else []
        if data:
            raise ValueError("Unexpected field")

        # Set / General
        (title, xmin, xmax, xlabel, xscale, ymin, ymax, ylabel, yscale, show_legend) = general

        if axes.get_xscale() != xscale:
            axes.set_xscale(xscale)
        if axes.get_yscale() != yscale:
            axes.set_yscale(yscale)

        axes.set_title(title)
        axes.set_xlim(xmin, xmax)
        axes.set_xlabel(xlabel)
        axes.set_ylim(ymin, ymax)
        axes.set_ylabel(ylabel)

        # Restore the unit data
        axes.xaxis.converter = xconverter
        axes.yaxis.converter = yconverter
        axes.xaxis.set_units(xunits)
        axes.yaxis.set_units(yunits)
        axes.xaxis._update_axisinfo()
        axes.yaxis._update_axisinfo()

        # Set / Curves
        for index, curve in enumerate(curves):
            line = labeled_lines[index][1]
            (
                label,
                linestyle,
                drawstyle,
                linewidth,
                color,
                marker,
                markersize,
                markerfacecolor,
                markeredgecolor,
            ) = curve
            line.set_label(label)
            line.set_linestyle(linestyle)
            line.set_drawstyle(drawstyle)
            line.set_linewidth(linewidth)
            rgba = mcolors.to_rgba(color)
            line.set_alpha(None)
            line.set_color(rgba)
            if marker != "none":
                line.set_marker(marker)
                line.set_markersize(markersize)
                line.set_markerfacecolor(markerfacecolor)
                line.set_markeredgecolor(markeredgecolor)

        # Set ScalarMappables.
        for index, mappable_settings in enumerate(mappables):
            mappable = labeled_mappables[index][1]
            if len(mappable_settings) == 5:
                label, cmap, low, high, interpolation = mappable_settings
                mappable.set_interpolation(interpolation)
            elif len(mappable_settings) == 4:
                label, cmap, low, high = mappable_settings
            mappable.set_label(label)
            mappable.set_cmap(cm.get_cmap(cmap))
            mappable.set_clim(*sorted([low, high]))

        # Set / Legend
        if show_legend:
            (
                leg_draggable,
                leg_ncol,
                leg_fontsize,
                leg_frameon,
                leg_framecolor,
                leg_bgcolor,
                leg_framealpha,
            ) = legend

            new_legend = axes.legend(
                ncol=leg_ncol,
                fontsize=float(leg_fontsize),
                frameon=leg_frameon,
                framealpha=leg_framealpha,
            )
            frame = new_legend.get_frame()
            frame.set_edgecolor(leg_framecolor)
            frame.set_facecolor(leg_bgcolor)
            new_legend.set_draggable(leg_draggable)
        else:
            axes.legend().set_visible(False)

        # Redraw
        figure = axes.get_figure()
        figure.canvas.draw()
        if not (axes.get_xlim() == orig_xlim and axes.get_ylim() == orig_ylim):
            figure.canvas.toolbar.push_current()

    _formlayout.fedit(
        datalist,
        title="Figure options",
        parent=parent,
        icon=QtGui.QIcon(str(cbook._get_data_path("images", "qt4_editor_options.svg"))),
        apply=apply_callback,
    )
