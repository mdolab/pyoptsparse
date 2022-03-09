# External modules
import Tkinter as Tk
import dill
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

matplotlib.use("TkAgg")

root = Tk.Tk()

# Load figure from disk and display
fig = dill.load(open("saved_figure.pickle", "rb"))

"""
The above code loads in the figure that was saved in OptView.
fig is a matplotlib object that can be altered and saved like any
regular figure.
The code at the bottom renders the image for immediate display.
Add your specific plot formatting code as necessary below this comment
string but before the bottom code.
"""

# Add customization code below
ax = fig.axes[0]
ax.set_title("Example title")
ax.set_ylabel("Example y-axis")


# Display the altered figure
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
toolbar = NavigationToolbar2TkAgg(canvas, root)
toolbar.update()
canvas.show()
Tk.mainloop()
