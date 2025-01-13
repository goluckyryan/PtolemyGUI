#!/usr/bin/python3

from PyQt6.QtWidgets import (
  QGridLayout, QWidget, QCheckBox
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar

from ExtractXsecPy import read_DWBA

# Set backend to a Qt-compatible one
plt.switch_backend('QtAgg')  # Or use 'Qt5Agg' if there are still issues

class MatPlotLibWindow(QWidget):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("DWBA Plot")
    self.resize(800, 600)

    self.x = []
    self.data = []
    self.headers = []

    self.figure, self.ax = plt.subplots()
    self.canvas = FigureCanvas(self.figure)
    self.toolbar = NavigationToolbar(self.canvas, self)

    self.showMarker_checkBox = QCheckBox("Show Markers")
    self.showMarker_checkBox.stateChanged.connect(self.plot_matplotlib_graph)

    self.log_scale_checkbox = QCheckBox("Use Log Scale for Y-Axis")
    self.log_scale_checkbox.setChecked(True)
    self.log_scale_checkbox.stateChanged.connect(self.plot_matplotlib_graph)

    self.gridline_checkbox = QCheckBox("Show Gridlines")
    self.gridline_checkbox.setChecked(True)
    self.gridline_checkbox.stateChanged.connect(self.plot_matplotlib_graph)


    layout = QGridLayout(self)
    layout.addWidget(self.toolbar, 0, 0, 1, 3)
    layout.addWidget(self.showMarker_checkBox, 1, 0)
    layout.addWidget(self.log_scale_checkbox, 1, 1)
    layout.addWidget(self.gridline_checkbox, 1, 2)
    layout.addWidget(self.canvas, 2, 0, 5, 3)

    self.ylabel = 'd.s.c.[mb/sr]'

  def set_ylable(self, newY_Label):
    self.ylabel = newY_Label

  def read_data(self,file_path):
    self.headers, self.x, self.data = read_DWBA(file_path)

  def plot_matplotlib_graph(self):
    self.ax.clear()

    plotStyle = '-' if not self.showMarker_checkBox.isChecked() else '-o'

    for i, y in enumerate(self.data):
      self.ax.plot(self.x, y, plotStyle, label=self.headers[i + 1])

    self.ax.set_xlabel(r"$\theta_{cm}$ [Deg]")
    self.ax.set_ylabel(self.ylabel)
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
    # self.figure.tight_layout()
    self.figure.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)


    self.canvas.draw_idle()

