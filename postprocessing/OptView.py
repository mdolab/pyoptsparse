"""

Provides interactive visualization of optimization results created by
pyOptSparse. Figures produced here can be saved as images or pickled
for future customization.

John Jasa 2015

"""

# ======================================================================
# Standard Python modules
# ======================================================================
import os
import argparse
import shelve
import tkFont
import Tkinter as Tk

# ======================================================================
# External Python modules
# ======================================================================
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,\
    NavigationToolbar2TkAgg
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import numpy
from pyoptsparse.pyoptsparse.sqlitedict.sqlitedict import SqliteDict

class Display(object):

    """
    Container for display parameters, properties, and objects.
    This includes a canvas for MPL plots and a bottom area with widgets.

    """

    def __init__(self, histFileName, outputDir):

        self.root = Tk.Tk()
        self.root.wm_title("OptView")

        try:
            icon_dir = os.path.dirname(os.path.abspath(__file__))
            icon_name = 'OptViewIcon.gif'
            icon_dir_full = os.path.join(icon_dir, icon_name)
            img = Tk.PhotoImage(file=icon_dir_full)
            self.root.tk.call('wm', 'iconphoto', self.root._w, img)
        except: # bare except because error is not in standard Python
            pass

        self.f = plt.figure(dpi=100, facecolor='white')

        self.canvas = FigureCanvasTkAgg(self.f, master=self.root)
        self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(self.canvas, self.root)
        toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        matplotlib.rcParams.update({'font.size': 16})

        self.arr_active = 0
        self.plots = []
        self.annotate = None
        self.histFileName = histFileName
        self.outputDir = outputDir

        self.OptimizationHistory()

    def OptimizationHistory(self):
        """
        Reads in database history file and stores contents.
        Function information is stored as a dict in func_data,
        variable information is stored as a dict in var_data,
        and bounds information is stored as a dict in bounds.
        """
        self.func_data = {}
        self.var_data = {}

        self.num_iter = 0

        try:
            db = shelve.open(self.histFileName, 'r')
        except: # bare except because error is not in standard Python
            db = SqliteDict(self.histFileName)
            
        nkey = int(db['last'])
        self.iter_type = numpy.zeros(nkey)

        # Check to see if there is bounds information in the hst file
        try:
            self.bounds = dict(
                db['varBounds'].items() + db['conBounds'].items())
        except KeyError:
            pass

        for i in xrange(nkey):
            key = '%d' % i
            keyp1 = '%d' % (i + 1)

            try:

                f = db[key]['funcs']
                try:
                    db[keyp1]['funcsSens']
                    self.iter_type[i] = 2
                except KeyError:
                    pass
                try:
                    db[keyp1]['funcs']
                    self.iter_type[i] = 1
                except KeyError:
                    pass

                for key in sorted(f):
                    if key not in self.func_data:
                        self.func_data[key] = []
                    if numpy.isscalar(f[key]):
                        self.func_data[key].append(f[key])
                    try:
                        if f[key].shape[0] > 1:
                            self.func_data[key].append(f[key])
                    except (IndexError, AttributeError):
                        pass

                try:
                    db[key]['funcsSens']
                except KeyError:
                    pass
                self.num_iter += 1

            except KeyError:
                pass

        for i in xrange(nkey):
            key = '%d' % i
            keyp1 = '%d' % (i + 1)
            if self.iter_type[i]:
                f = db[key]['xuser']
                for key in sorted(f):
                    if key not in self.var_data:
                        self.var_data[key] = []
                    if numpy.isscalar(f[key]):
                        self.var_data[key].append(f[key])
                    try:
                        if f[key].shape[0] > 1:
                            self.var_data[key].append(f[key])
                    except IndexError:
                        pass

    def quit(self):
        """
        Destroy GUI window cleanly if quit button pressed.
        """
        self.root.quit()
        self.root.destroy()

    def error_display(self, string="That option not supported"):
        """
        Display error string on canvas when invalid options selected.
        """
        self.f.clf()
        a = self.f.add_subplot(111)
        a.text(
            0.05,
            .9,
            "Error: " + string,
            fontsize=20,
            transform=a.transAxes)
        self.canvas.show()

    def warning_display(self, string="That option not supported"):
        """
        Display warning message on canvas as necessary.
        """
        a = plt.gca()
        a.text(
            0.05,
            1.04,
            "Warning: " + string,
            fontsize=20,
            transform=a.transAxes)
        self.canvas.show()

    def plot_bounds(self, val, a, color):
        """
        Plot the bounds if selected.
        """

        if val not in self.bounds:
            for ii, char in enumerate(reversed(val)):
                if char == '_':
                    split_loc = len(val) - ii
                    break
            val_name = val[:split_loc - 1]
            val_num = int(val[split_loc:])
            lower = [self.bounds[val_name]['lower'][val_num]]
            upper = [self.bounds[val_name]['upper'][val_num]]
        else:
            lower = self.bounds[val]['lower']
            upper = self.bounds[val]['upper']
            if disp.var.get() == 0:
                color = None

        lower = list(lower)
        upper = list(upper)

        a.margins(None, .02)
        a.set_color_cycle(color)
        for lower_bound in lower:
            if lower_bound is not None and abs(lower_bound) < 1e18:
                a.plot(
                    [0, self.num_iter - 1], [
                        lower_bound, lower_bound],
                    "--", linewidth=2
                )

        a.set_color_cycle(color)
        for upper_bound in upper:
            if upper_bound is not None and abs(upper_bound) < 1e18:
                a.plot(
                    [0, self.num_iter - 1], [
                        upper_bound, upper_bound],
                    "--", label=val + ' bounds', linewidth=2)

    def orig_plot(self, dat, val, values, a, i=0):
        """
        Plots the original data values from the history file.
        """
        color = None
        try:
            array_size = len(dat[val][0])
            if self.var_minmax.get():
                minmax_list = []
                for minmax in dat[val]:
                    minmax_list.append(
                        [numpy.min(minmax), numpy.max(minmax)])
                plots = a.plot(
                    minmax_list,
                    "o-",
                    label=val,
                    markeredgecolor='none')

            elif array_size < 20 or self.var_showall.get():
                plots = a.plot(
                    dat[val],
                    "o-",
                    label=val,
                    markeredgecolor='none')

                a.set_ylabel(val)
                self.color_error_flag = 1
            else:
                self.error_display("Too many values to display")

        except TypeError:
            if self.var.get() == 0:
                try:
                    if self.var_bounds.get() and val not in self.bounds:
                        cc = (
                            matplotlib.rcParams['axes.color_cycle'] * 10
                        )
                        color = cc[i]
                        a.set_color_cycle(color)
                except UnboundLocalError:
                    pass
            else:
                a.set_ylabel(val)
            plots = a.plot(
                dat[val],
                "o-",
                label=val,
                markeredgecolor='none')
        try:
            if len(plots) > 1:
                for i, plot in enumerate(plots):
                    self.plots.append([plot, i])
            else:
                self.plots.append([plots[0], -1])
        except UnboundLocalError:
            self.error_display("Too many values to display")
        try:
            if self.var_bounds.get():
                self.plot_bounds(val, a, color)

        except (UnboundLocalError, ValueError):
            if len(values) > 1:
                pass
            else:
                self.error_display("No bounds information")

    def plot_selected(self, values, dat):
        """
        Plot data based on selected keys within listboxes.
        """
        fail = 0
        self.color_error_flag = 0
        self.f.clf()
        self.plots = []
        try:
            if self.var_bounds.get():
                try:
                    self.bounds
                except AttributeError:
                    self.error_display("No bounds information in history file")
                    fail = 1

            # Plot on shared axes
            if self.var.get() == 0 and not fail:
                a = self.f.add_subplot(111)

                # Calculate and plot the delta values if selected
                if self.var_del.get():
                    for idx, val in enumerate(values):
                        newdat = []
                        for i, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][i - 2]))
                        plots = a.plot(
                            range(1, self.num_iter),
                            newdat[1:],
                            "o-",
                            label=val,
                            markeredgecolor='none')
                        if len(plots) > 1:
                            for i, plot in enumerate(plots):
                                self.plots.append([plot, i])
                        else:
                            self.plots.append([plots[0], -1])

                # Otherwise plot original data
                else:
                    for i, val in enumerate(values):
                        self.orig_plot(dat, val, values, a, i)

                if i > 0 and self.color_error_flag and self.var_bounds.get():
                    self.warning_display(
                        "Line color for bounds may not match data color")

                # Plot using log scale if selected
                if self.var_log.get():
                    a.set_yscale('log')

                if self.var_legend.get():
                    a.legend(loc='best')

                plt.subplots_adjust(right=.95)
                a.set_xlabel('iteration')
                a.set_xlim(0, self.num_iter - 1)
                self.canvas.show()

            # Plot on individual vertical axes
            elif self.var.get() == 1 and not fail:

                # Set window sizing parameters for when additional axes are
                # added
                n = len(values)
                plt.figure(self.f.number)
                par_list = [[] for i in range(n)]  # make into array
                par_list[0] = host_subplot(111, axes_class=AA.Axes)
                size_list = [.95, .95, .93, .83, .73, .63, .53, .43, .33]
                plt.subplots_adjust(right=size_list[n])

                for i in range(1, n):
                    par_list[i] = par_list[0].twinx()

                offset = 60
                for i in range(2, n):
                    new_fixed_axis = par_list[
                        i].get_grid_helper().new_fixed_axis
                    par_list[i].axis["right"] = new_fixed_axis(
                        loc="right", axes=par_list[i], offset=(offset * i ** 1.15, 0))
                    par_list[i].axis["right"].toggle(all=True)

                p_list = [[] for i in range(n)]

                # Compute and plot delta values if selected
                if self.var_del.get():
                    for i, val in enumerate(values):
                        newdat = []
                        for idx, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][idx - 2]))
                        p_list[i], = par_list[i].plot(range(1, self.num_iter),
                                                    newdat[1:], "o-", label=val,
                                                    markeredgecolor='none')
                        par_list[i].set_ylabel(val)
                # Otherwise plot original data
                else:
                    for i, val in enumerate(values):
                        cc = (matplotlib.rcParams['axes.color_cycle'] * 10)
                        par_list[i].set_color_cycle(cc[i])
                        p_list[i], = par_list[i].plot(
                            dat[val], "o-", label=val, markeredgecolor='none')
                        par_list[i].set_ylabel(val)

                        try:
                            if self.var_bounds.get():
                                self.plot_bounds(val, par_list[i], cc[i])
                        except (ValueError, UnboundLocalError):
                            if len(values) > 1:
                                pass
                            else:
                                self.error_display("No bounds information")

                # Plot using log scale if selected
                if self.var_log.get():
                    for ax in par_list:
                        ax.set_yscale('log')

                par_list[0].set_xlim(0, self.num_iter - 1)
                par_list[0].set_xlabel('iteration')
                if self.var_legend.get():
                    par_list[0].legend(loc='best')
                for i, plot in enumerate(p_list):
                    self.plots.append([plot, i])

                self.canvas.show()

            # Plot on stacked axes with shared x-axis
            elif self.var.get() == 2 and not fail:
                n = len(values)

                # Compute and plot delta values if selected
                if self.var_del.get():
                    a = []
                    for i, val in enumerate(values):
                        newdat = []
                        for idx, value in enumerate(dat[val], start=1):
                            newdat.append(abs(value - dat[val][idx - 2]))
                        a.append(self.f.add_subplot(n, 1, i + 1))
                        plots = a[i].plot(range(1, self.num_iter), newdat[1:],
                                          "o-", label=val, markeredgecolor='none')
                        a[i].set_ylabel('delta ' + val)
                        self.plots.append([plots[0], -1])

                # Otherwise plot original data
                else:
                    a = []
                    for i, val in enumerate(values):
                        a.append(self.f.add_subplot(n, 1, i + 1))
                        self.orig_plot(dat, val, values, a[i])

                # Plot using log scale if selected
                if self.var_log.get():
                    for ax in a:
                        ax.set_yscale('log')

                # Turn off horiztonal axes above the bottom plot
                a[-1].set_xlabel('iteration')
                for ax in a:
                    if ax != a[-1]:
                        ax.spines['bottom'].set_visible(False)
                        ax.set_xticklabels([])
                        ax.xaxis.set_major_locator(plt.NullLocator())
                    ax.spines['top'].set_visible(False)
                    for tic in ax.xaxis.get_major_ticks():
                        tic.tick2On = False
                    ax.tick_params(
                        axis='y',
                        which='both',
                        labelleft='off',
                        labelright='on')
                    ax.set_xlim(0, self.num_iter - 1)

                plt.subplots_adjust(right=.95)
                self.canvas.show()
        except ValueError:
            self.error_display()

    def onselect(self, evt, data_name):
        """
        Update current plot with selected data from listboxes.
        Also checks if the data is array-type and provides an
        additional listbox to select data within that array.
        """
        w = evt.widget
        values = [w.get(int(i)) for i in w.curselection()]
        self.plot_selected(values, data_name)
        if len(values) == 1:
            try:
                data = data_name[values[0]]
                data[0][0]
                self.v.set(values[0])
                self.lb_arr.delete(0, Tk.END)
                for i, val in enumerate(data[0]):
                    self.lb_arr.insert(Tk.END, values[0] + '_' + str(i))
                self.arr_title.pack(side=Tk.TOP)
                self.scrollbar_arr.pack(side=Tk.RIGHT, fill=Tk.Y)
                self.lb_arr.pack(side=Tk.RIGHT)
                self.arr_active = 1
            except IndexError:
                self.lb_arr.pack_forget()
                self.scrollbar_arr.pack_forget()
                self.arr_title.pack_forget()
                self.arr_active = 0

    def onselect_arr(self, evt):
        """
        Obtain the selected plotting values from the array-based variable listbox.
        """
        w = evt.widget
        values = [int(i) for i in w.curselection()]

        # Get the currently selected functions/variables
        func_sel = self.lb_func.curselection()
        var_sel = self.lb_var.curselection()
        if len(func_sel):
            values_orig = [self.lb_func.get(i) for i in func_sel]
            dat = self.func_data[values_orig[0]]
        elif len(var_sel):
            values_orig = [self.lb_var.get(i) for i in var_sel]
            dat = self.var_data[values_orig[0]]

        # Add the array-based information to the listbox for selection
        self.arr_data = {}
        self.val_names = []
        for i, val in enumerate(values):
            self.val_names.append(values_orig[0] + '_{0}'.format(val))
            self.arr_data[self.val_names[i]] = []
            for ind_dat in dat:
                self.arr_data[self.val_names[i]].append(ind_dat[val])
        self.plot_selected(self.val_names, self.arr_data)

    def update_graph(self):
        """
        Produce an updated graph based on user options.
        """
        if self.var_minmax.get() and self.var_showall.get():
            self.error_display("Cannot show all and min/max at same time")
        else:
            func_sel = self.lb_func.curselection()
            var_sel = self.lb_var.curselection()
            arr_sel = self.lb_arr.curselection()
            if len(arr_sel) and self.arr_active:
                self.plot_selected(self.val_names, self.arr_data)
            elif len(func_sel):
                values = [self.lb_func.get(i) for i in func_sel]
                self.plot_selected(values, self.func_data)
            elif len(var_sel):
                values = [self.lb_var.get(i) for i in var_sel]
                self.plot_selected(values, self.var_data)

    def save_figure(self):
        """
        Save the current figure using the selected variables as the filename.
        """
        func_sel = self.lb_func.curselection()
        var_sel = self.lb_var.curselection()
        arr_sel = self.lb_arr.curselection()
        values = []
        if len(arr_sel) and self.arr_active:
            values = self.val_names
        elif len(func_sel):
            values = [self.lb_func.get(i) for i in func_sel]
        elif len(var_sel):
            values = [self.lb_var.get(i) for i in var_sel]
        groups = ''
        for string in values:
            groups += string + '_'
        fname = groups + '.png'
        fpathname = os.path.join(self.outputDir, fname)
        plt.savefig(fpathname)
        fname = 'saved_figure.pickle'
        fpathname = os.path.join(self.outputDir, fname)
        try:
            import dill
            dill.dump(self.f, file(fpathname, 'wb'))
        except ImportError:
            pass            

    def save_all_figues(self):
        """
        Batch save all individual figures from functions and variables.
        """
        for data_name in [self.func_data, self.var_data]:
            for key in data_name:
                fig = plt.figure()
                plt.plot(data_name[key], 'ko-')
                plt.title(key)
                fname = key + '.png'
                fpathname = os.path.join(self.outputDir, fname)
                plt.savefig(fpathname)
                plt.clf()

    def var_search(self, _):
        """
        Remove listbox entries that do not contain user-inputted string,
        used to search through outputted data.
        """
        self.lb_func.delete(0, Tk.END)
        self.lb_var.delete(0, Tk.END)
        for key in sorted(self.func_data):
            self.lb_func.insert(Tk.END, key)
        for key in sorted(self.var_data):
            self.lb_var.insert(Tk.END, key)

        search_entry = self.entry_search.get()
        func_range = range(len(self.func_data))
        for i in func_range[::-1]:
            if search_entry.lower() not in self.lb_func.get(i).lower():
                self.lb_func.delete(i)

        var_range = range(len(self.var_data))
        for i in var_range[::-1]:
            if search_entry.lower() not in self.lb_var.get(i).lower():
                self.lb_var.delete(i)

        if not self.lb_var.get(1) and not self.lb_func.get(1):
            if self.lb_var.get(0):
                self.lb_var.select_set(0)
            if self.lb_func.get(0):
                self.lb_func.select_set(0)
            self.update_graph()

    def update_font(self, val):
        """
        Set the font for matplotlib based on slider.
        """
        matplotlib.rcParams.update({'font.size': int(val)})
        self.update_graph()

    def refresh_history(self):
        """
        Refresh opt_his data if the history file has been updated.
        """
        old_funcs = []
        for key in self.func_data:
            old_funcs.append(key)
        old_vars = []
        for key in self.var_data:
            old_vars.append(key)

        self.OptimizationHistory()
        self.update_graph()

        new_funcs = []
        for key in self.func_data:
            new_funcs.append(key)
        new_vars = []
        for key in self.var_data:
            new_vars.append(key)

        if not (old_funcs == new_funcs and old_vars == new_vars):
            self.var_search('dummy')

    def update_all(self):
        """
        Updates the history and graph.
        """
        self.refresh_history(self.plotAll.get())
        self.update_graph()

    def on_move(self, event):
        """
        Checks to see if the cursor is over a plot and provides a
        hovering label if necessary.
        """
        try:
            self.annotation.remove()
        except (AttributeError, ValueError):
            pass
        visibility_changed = False
        point_selected = None
        for point in self.plots:
            if point[0].contains(event)[0]:
                point_selected = point

        if point_selected:
            visibility_changed = True
            ax = point_selected[0].get_axes()
            label = point_selected[0].get_label()
            if point_selected[1] >= 0:
                label = label + '_' + str(point_selected[1])
            self.annotation = ax.annotate(label,
                                          xy=(event.xdata,
                                              event.ydata), xycoords='data',
                                          xytext=(
                                          event.xdata, event.ydata), textcoords='data',
                                          horizontalalignment="left",
                                          bbox=dict(
                                          boxstyle="round", facecolor="w",
                                          edgecolor="0.5", alpha=0.8),
                                          )
        else:
            try:
                self.annotation.remove()
            except (AttributeError, ValueError):
                pass

        self.canvas.show()

    def draw_GUI(self):
        """
        Create the frames and widgets in the bottom section of the canvas.
        """
        font = tkFont.Font(family="Helvetica", size=10)

        sel_frame = Tk.Frame(self.root)
        sel_frame.pack(side=Tk.LEFT)

        # Produce a frame and listbox to contain function information
        func_frame = Tk.Frame(sel_frame)
        func_frame.pack(side=Tk.LEFT, fill=Tk.Y, padx=20)
        func_title = Tk.Label(func_frame, text="Functions", font=font)
        func_title.pack(side=Tk.TOP)
        self.lb_func = Tk.Listbox(
            func_frame,
            name='lb_func',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30)
        self.lb_func.pack(side=Tk.LEFT)
        for key in sorted(self.func_data):
            self.lb_func.insert(Tk.END, key)
        scrollbar_func = Tk.Scrollbar(func_frame)
        scrollbar_func.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lb_func.config(yscrollcommand=scrollbar_func.set)
        scrollbar_func.config(command=self.lb_func.yview)
        self.lb_func.bind(
            '<<ListboxSelect>>',
            lambda event: self.onselect(
                event,
                self.func_data))

        # Produce a frame and listbox to contain variable information
        var_frame = Tk.Frame(sel_frame)
        var_frame.pack(side=Tk.RIGHT, fill=Tk.Y, padx=20)
        var_title = Tk.Label(var_frame, text="Variables", font=font)
        var_title.pack(side=Tk.TOP)
        scrollbar_var = Tk.Scrollbar(var_frame)
        scrollbar_var.pack(side=Tk.RIGHT, fill=Tk.Y)
        self.lb_var = Tk.Listbox(
            var_frame,
            name='lb_var',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30)
        self.lb_var.pack(side=Tk.RIGHT)
        for key in sorted(self.var_data):
            self.lb_var.insert(Tk.END, key)
        self.lb_var.config(yscrollcommand=scrollbar_var.set)
        scrollbar_var.config(command=self.lb_var.yview)
        self.lb_var.bind(
            '<<ListboxSelect>>',
            lambda event: self.onselect(
                event,
                self.var_data))

        # Produce a frame and listbox to contain array-based variable
        # information
        arr_frame = Tk.Frame(self.root)
        arr_frame.pack(side=Tk.RIGHT, fill=Tk.Y, padx=20)
        self.v = Tk.StringVar()
        self.arr_title = Tk.Label(
            arr_frame,
            text="Array variables",
            font=font,
            textvariable=self.v)
        self.scrollbar_arr = Tk.Scrollbar(arr_frame)
        self.lb_arr = Tk.Listbox(
            arr_frame,
            name='lb_arr',
            selectmode=Tk.EXTENDED,
            font=font,
            width=30,
            exportselection=0)
        self.lb_arr.config(yscrollcommand=self.scrollbar_arr.set)
        self.scrollbar_arr.config(command=self.lb_arr.yview)
        self.lb_arr.bind('<<ListboxSelect>>', self.onselect_arr)

        # Create options frame with buttons and checkboxes for user
        options_frame = Tk.Frame(self.root)
        options_frame.pack(side=Tk.LEFT, fill=Tk.Y, pady=10)

        button1 = Tk.Button(
            options_frame,
            text='Refresh history',
            command=self.refresh_history,
            font=font)
        button1.grid(row=0, column=2, padx=5, sticky=Tk.W)

        button2 = Tk.Button(
            options_frame,
            text='Save all figures',
            command=self.save_all_figues,
            font=font)
        button2.grid(row=1, column=2, padx=5, sticky=Tk.W)

        button3 = Tk.Button(
            options_frame,
            text='Save figure',
            command=self.save_figure,
            font=font)
        button3.grid(row=2, column=2, padx=5, sticky=Tk.W)

        button4 = Tk.Button(
            options_frame,
            text='Quit',
            command=self.quit,
            font=font)
        button4.grid(row=3, column=2, padx=5, sticky=Tk.W)

        # Plot options
        self.var = Tk.IntVar()
        c1 = Tk.Radiobutton(
            options_frame, text="Shared axes", variable=self.var,
            command=self.update_graph, font=font, value=0)
        c1.grid(row=0, column=0, sticky=Tk.W)

        c2 = Tk.Radiobutton(
            options_frame, text="Multiple axes", variable=self.var,
            command=self.update_graph, font=font, value=1)
        c2.grid(row=1, column=0, sticky=Tk.W)

        c3 = Tk.Radiobutton(
            options_frame, text="Stacked plots", variable=self.var,
            command=self.update_graph, font=font, value=2)
        c3.grid(row=2, column=0, sticky=Tk.W)

        self.var_del = Tk.IntVar()
        c4 = Tk.Checkbutton(
            options_frame,
            text="Absolute delta values",
            variable=self.var_del,
            command=self.update_graph,
            font=font)
        c4.grid(row=0, column=1, sticky=Tk.W)

        self.var_log = Tk.IntVar()
        c5 = Tk.Checkbutton(
            options_frame,
            text="Log scale",
            variable=self.var_log,
            command=self.update_graph,
            font=font)
        c5.grid(row=1, column=1, sticky=Tk.W)

        # Option to only show the min and max of array-type variables
        self.var_minmax = Tk.IntVar()
        c6 = Tk.Checkbutton(
            options_frame,
            text="Min/max for arrays",
            variable=self.var_minmax,
            command=self.update_graph,
            font=font)
        c6.grid(row=2, column=1, sticky=Tk.W)

        # Option to show all values for array-type variables
        self.var_showall = Tk.IntVar()
        c7 = Tk.Checkbutton(
            options_frame,
            text="Show all for arrays",
            variable=self.var_showall,
            command=self.update_graph,
            font=font)
        c7.grid(row=3, column=1, sticky=Tk.W)

            # Option to show legend
        self.var_legend = Tk.IntVar()
        c8 = Tk.Checkbutton(
            options_frame,
            text="Show legend",
            variable=self.var_legend,
            command=self.update_graph,
            font=font)
        self.var_legend.set(1)
        c8.grid(row=4, column=1, sticky=Tk.W)

        # Option to show bounds
        self.var_bounds = Tk.IntVar()
        c9 = Tk.Checkbutton(
            options_frame,
            text="Show bounds",
            variable=self.var_bounds,
            command=self.update_graph,
            font=font)
        c9.grid(row=5, column=1, sticky=Tk.W)

        lab = Tk.Label(
            options_frame,
            text="Search for a function/variable:",
            font=font)
        lab.grid(row=6, column=0, columnspan=2, pady=10, sticky=Tk.W)

        # Search box to filter displayed functions/variables
        vs = Tk.StringVar()
        vs.trace("w", lambda name, index, mode, vs=vs: self.var_search(vs))
        self.entry_search = Tk.Entry(
            options_frame, text="Search", textvariable=vs,
            font=font)
        self.entry_search.grid(row=6, column=2, pady=10, sticky=Tk.W)

        lab_font = Tk.Label(
            options_frame,
            text="Font size for plots:",
            font=font)
        lab_font.grid(row=7, column=0, sticky=Tk.S)

        w = Tk.Scale(options_frame, from_=6, to=24, orient=Tk.HORIZONTAL,
                     resolution=2, command=self.update_font, font=font)
        w.set(16)
        w.grid(row=7, column=1)

if __name__ == '__main__':
    # Called only if this script is run as main.

    # ======================================================================
    # Input Information
    # ======================================================================
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'histFile', nargs='?', type=str, default='opt_hist.hst',
        help="Specify the history file to be plotted")
    parser.add_argument('output', nargs='?', type=str, default='./',
                        help="Specify the output directory")
    args = parser.parse_args()
    histFileName = args.histFile
    outputDir = args.output

    # Check that the output directory is available. Create it if not
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
    # Initialize display parameters, obtain history, and draw GUI
    disp = Display(histFileName, outputDir)
    disp.draw_GUI()
    disp.root.protocol("WM_DELETE_WINDOW", disp.quit)
    on_move_id = disp.f.canvas.mpl_connect('motion_notify_event', disp.on_move)
    Tk.mainloop()
