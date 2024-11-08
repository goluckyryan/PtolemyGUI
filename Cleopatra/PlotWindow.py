#!/usr/bin/python3

from PyQt6.QtWidgets import (
  QGridLayout, QWidget, QCheckBox
)

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from ExtractXsecPy import read_DWBA

class PlotWindow(QWidget):
  def __init__(self, windowTitle):
    super().__init__()

    self.setWindowTitle(windowTitle)
    self.resize(800, 600)

    self.default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    self.log_scale_checkbox = QCheckBox("Use Log Scale for Y-Axis")
    self.log_scale_checkbox.setChecked(True)
    self.log_scale_checkbox.stateChanged.connect(self.plot_graph)

    self.gridline_checkbox = QCheckBox("Show Gridlines")
    self.gridline_checkbox.stateChanged.connect(self.plot_graph)

    self.showMarker_checkBox = QCheckBox("Show Markers")
    self.showMarker_checkBox.stateChanged.connect(self.plot_graph)

    self.figure, self.ax = plt.subplots()
    self.canvas = FigureCanvas(self.figure)
    self.toolbar = NavigationToolbar(self.canvas, self)

    layout = QGridLayout(self)
    layout.addWidget(self.toolbar, 0, 0, 1, 3)
    layout.addWidget(self.showMarker_checkBox, 1, 0)
    layout.addWidget(self.log_scale_checkbox, 1, 1)
    layout.addWidget(self.gridline_checkbox, 1, 2)
    layout.addWidget(self.canvas, 2, 0, 5, 3)

    self.xData = []
    self.yData_list = []
    self.header_list = []
    self.yTitle = ""

    self.x_exp = [] # x positions
    self.x_err = [] # x uncertainties (errors)
    self.y_exp = [] # y positions
    self.y_err = [] # y errors
    self.dataName = ""
    self.fitOption = []

    self.para = [] # fit parameters
    self.perr = [] # fit error
    self.chi_square = [] # fit Chi-squared

  def set_plot_data(self, xData, yData_list, header_list, yTitle):
    self.xData = xData
    self.yData_list = yData_list
    self.header_list = header_list
    self.yTitle = yTitle

  def set_expData(self, expData, fitOption, dataName_list, ID):
    self.x_exp = expData[ID][:, 0]
    self.x_err = expData[ID][:, 1]
    self.y_exp = expData[ID][:, 2]
    self.y_err = expData[ID][:, 3]
    self.dataName = dataName_list[ID]
    self.fitOption = fitOption[ID]

  def read_Xsec(self, file_path):
    headers, dataX, data = read_DWBA(file_path)
    self.xData = dataX
    self.yData_list = data
    self.header_list = headers[1:]

  def set_fitResult(self, para, perr, chi_sq):
    self.para = para
    self.perr = perr
    self.chi_square = chi_sq

  def plot_graph(self):
    self.ax.clear()

    plotStyle = '-' if not self.showMarker_checkBox.isChecked() else '-o'

    for i, y in enumerate(self.yData_list):
      self.ax.plot(self.xData, y, plotStyle, label=self.header_list[i])

    self.ax.set_xlabel('Angle_CM [deg]')
    self.ax.set_ylabel(self.yTitle)
    self.ax.legend(loc='upper right', frameon=True)

    # Apply log scale for y-axis if selected
    if self.log_scale_checkbox.isChecked():
      self.ax.set_yscale('log')
    else:
      self.ax.set_yscale('linear')

    self.ax.autoscale(enable=True, axis='x', tight=True)
    self.figure.tight_layout()

  def plot_Fit(self):
    self.ax.clear()

    self.ax.errorbar(self.x_exp, self.y_exp, xerr=self.x_err, yerr=self.y_err, 
                    fmt='x', label='Experimental Data', color='black', markersize = 15, elinewidth=2)

    self.ax.set_xlabel('Angle_CM [deg]')
    self.ax.set_ylabel(self.yTitle)
    self.ax.legend(loc='upper right', frameon=True)

    # Apply log scale for y-axis if selected
    if self.log_scale_checkbox.isChecked():
      self.ax.set_yscale('log')
    else:
      self.ax.set_yscale('linear')

    self.ax.autoscale(enable=True, axis='x', tight=True)
    self.figure.tight_layout()

    for k in range(len(self.fitOption)):
      fitTheory = []
      fitTheory_lower = []
      fitTheory_upper = []

      xsecIDStr = self.fitOption[k].strip()
      xsecID = [int(part.strip()) for part in xsecIDStr.split('+')] if '+' in xsecIDStr else [int(xsecIDStr)]

      fitTheory.append(np.zeros_like(self.xData))
      for p, id in enumerate(xsecID):
        fitTheory += self.para[p] * np.interp(self.xData, self.xData, self.yData_list[id])

      fitTheory_upper.append(np.zeros_like(self.xData))
      fitTheory_lower.append(np.zeros_like(self.xData))
      
      for p, id in enumerate(xsecID):
        fitTheory_upper += (self.para[p] + self.perr[p]) * np.interp(self.xData, self.xData, self.yData_list[id])
        fitTheory_lower += (self.para[p] - self.perr[p]) * np.interp(self.xData, self.xData, self.yData_list[id])

    # Replace plt.title() with plt.text() to position the title inside the plot
    self.ax.text(0.05, 0.05, f'Fit for Exp Data : {self.dataName}', transform=plt.gca().transAxes,
        fontsize=12, verticalalignment='bottom', horizontalalignment='left', color='black')

    for i, fit in enumerate(fitTheory):
      self.ax.plot(self.xData, fit, label=f'Chi2:{self.chi_square[i]:.3f} | Xsec:{self.fitOption[i]}')
      self.ax.fill_between(self.xData, fitTheory_lower[i], fitTheory_upper[i], alpha=0.2)

    for i, _ in enumerate(self.para):
      self.ax.text(0.05, 0.1 + 0.05*i, f"Xsec-{self.fitOption[i].strip()}: {', '.join([f'{x:.3f}' for x in self.para[i]])} +/- {', '.join([f'{x:.3f}' for x in self.perr[i]])}" ,
                transform=plt.gca().transAxes, fontsize=12, 
                verticalalignment='bottom', horizontalalignment='left', color=self.default_colors[i])

    self.canvas.draw_idle()






