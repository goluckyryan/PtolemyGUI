#!/usr/bin/python3

from PyQt6.QtWidgets import (
  QGridLayout, QWidget, QCheckBox
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib import get_backend

# Set backend to a Qt-compatible one
plt.switch_backend('QtAgg')  # Or use 'Qt5Agg' if there are still issues

class MatPlotLibWindow(QWidget):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("DWBA Plot")
    self.setGeometry(100, 100, 800, 600)

    self.x = []
    self.data = []
    self.headers = []

    self.figure, self.ax = plt.subplots()
    self.canvas = FigureCanvas(self.figure)
    self.toolbar = NavigationToolbar(self.canvas, self)

    self.log_scale_checkbox = QCheckBox("Use Log Scale for Y-Axis")
    self.log_scale_checkbox.setChecked(True)
    self.log_scale_checkbox.stateChanged.connect(self.plot_matplotlib_graph)

    self.gridline_checkbox = QCheckBox("Show Gridlines")
    self.gridline_checkbox.stateChanged.connect(self.plot_matplotlib_graph)

    self.showMarker_checkBox = QCheckBox("Show Markers")
    self.showMarker_checkBox.setChecked(True)
    self.showMarker_checkBox.stateChanged.connect(self.plot_matplotlib_graph)

    layout = QGridLayout(self)
    layout.addWidget(self.toolbar, 0, 0, 1, 3)
    layout.addWidget(self.showMarker_checkBox, 1, 0)
    layout.addWidget(self.log_scale_checkbox, 1, 1)
    layout.addWidget(self.gridline_checkbox, 1, 2)
    layout.addWidget(self.canvas, 2, 0, 5, 3)

  def read_data(self,file_path):
    self.x = []  # List for the first column
    self.data = [] # 2D list for other columns
    self.headers = []  # List to store headers

    with open(file_path, 'r') as file:
      header_found = False  # Flag to indicate if the header has been found
      for line in file:
        # Skip lines that start with '#' and empty lines
        if line.startswith('#') or not line.strip():
          continue
        
        if not header_found:
          self.headers = line.split()  # Use the split parts as headers
          header_found = True  # Set the flag to True to skip this block in future iterations
          # print(f"ELab parts found: {elab_parts}")  # Print or process this as needed
          continue
        
        # Split the line by whitespace
        parts = line.split()
        if len(parts) > 0:  # Make sure there is at least one column
          self.x.append(float(parts[0]))  # First column
          # Append the rest of the columns to data
          if len(self.data) == 0:
            # Initialize the data array with the right number of sublists
            self.data = [[] for _ in range(len(parts) - 1)]
          for i in range(len(parts) - 1):
            self.data[i].append(float(parts[i + 1]))  # Rest of the columns
            
      # print(self.headers)

  def plot_matplotlib_graph(self):
    self.ax.clear()

    plotStyle = '-' if not self.showMarker_checkBox.isChecked() else '-o'

    for i, y in enumerate(self.data):
      self.ax.plot(self.x, y, plotStyle, label=self.headers[i + 1])

    self.ax.set_xlabel("Angle_CM [Deg]")
    self.ax.set_ylabel("Xsec [mb/sr]")
    self.ax.legend(loc='upper right', frameon=True)

    # Apply log scale for y-axis if selected
    if self.log_scale_checkbox.isChecked():
        self.ax.set_yscale('log')
    else:
        self.ax.set_yscale('linear')

    # Toggle gridlines
    if self.gridline_checkbox.isChecked():
      self.ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
    else:
      self.ax.grid(False)

    self.ax.autoscale(enable=True, axis='x', tight=True)
    self.figure.tight_layout()

    self.canvas.draw_idle()

