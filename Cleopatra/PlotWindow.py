#!/usr/bin/python3

from PyQt6.QtWidgets import (
  QGridLayout, QWidget, QCheckBox
)

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

# Set backend to a Qt-compatible one
plt.switch_backend('QtAgg')  # Or use 'Qt5Agg' if there are still issues

class FitPlotWindow(QWidget):
  def __init__(self, windowTitle):
    super().__init__()

    self.setWindowTitle(windowTitle)
    self.resize(800, 600)

    self.default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    self.log_scale_checkbox = QCheckBox("Use Log Scale for Y-Axis")
    self.log_scale_checkbox.setChecked(True)
    self.log_scale_checkbox.stateChanged.connect(self.plot_Fit)

    self.gridline_checkbox = QCheckBox("Show Gridlines")
    self.gridline_checkbox.stateChanged.connect(self.plot_Fit)

    self.legend_checkbox = QCheckBox("Show Legends")
    self.legend_checkbox.setChecked(True)
    self.legend_checkbox.stateChanged.connect(self.plot_Fit)

    self.figure, self.ax = plt.subplots()
    self.canvas = FigureCanvas(self.figure)
    self.toolbar = NavigationToolbar(self.canvas, self)

    layout = QGridLayout(self)
    layout.addWidget(self.toolbar, 0, 0, 1, 3)
    layout.addWidget(self.log_scale_checkbox, 1, 0)
    layout.addWidget(self.gridline_checkbox, 1, 1)
    layout.addWidget(self.legend_checkbox, 1, 2)
    layout.addWidget(self.canvas, 2, 0, 5, 3)

    self.setLayout(layout)

  def set_data(self, ID, expData, fitOption, dataName_list, xData, yData_list, headers, para, perr, chi_square ):
    self.x_exp = expData[ID][:, 0]
    self.x_err = expData[ID][:, 1]
    self.y_exp = expData[ID][:, 2]
    self.y_err = expData[ID][:, 3]
    self.dataName = dataName_list
    self.fitOption = fitOption[ID]

    self.para = para
    self.perr = perr
    self.chi_square = chi_square

    self.xData = xData
    self.yData_list = yData_list

    self.headers = headers

    print(self.dataName)
    print(self.fitOption)
    print(self.headers)

  def plot_Fit(self):
    self.ax.clear()

    fitTheory = []
    fitTheory_lower = []
    fitTheory_upper = []

    fitXsecID = []
    fitHeaders = []
    
    for k in range(len(self.fitOption)):

      xsecIDStr = self.fitOption[k].strip()
      xsecID = [int(part.strip()) for part in xsecIDStr.split('+')] if '+' in xsecIDStr else [int(xsecIDStr)]

      fitTheory.append(np.zeros_like(self.xData))
      fitTheory_upper.append(np.zeros_like(self.xData))
      fitTheory_lower.append(np.zeros_like(self.xData))
      
      for i, id in enumerate(xsecID):
        fitXsecID.append(id)
        fitHeaders.append(self.headers[id])
        fitTheory[k] += self.para[k][i] * np.interp(self.xData, self.xData, self.yData_list[id])
        fitTheory_upper[k] += (self.para[k][i] + self.perr[k][i]) * np.interp(self.xData, self.xData, self.yData_list[id])
        fitTheory_lower[k] += (self.para[k][i] - self.perr[k][i]) * np.interp(self.xData, self.xData, self.yData_list[id])


    for i, fit in enumerate(fitTheory):
      self.ax.plot(self.xData, fit, label=f'Chi2:{self.chi_square[i]:.3f} | Xsec:{self.fitOption[i]}')
      self.ax.fill_between(self.xData, fitTheory_lower[i], fitTheory_upper[i], alpha=0.2)

      if self.legend_checkbox.isChecked() :
        self.ax.text(0.98, 0.98 - 0.05*i, rf"Fit-{self.fitOption[i].strip()}: {', '.join([f'{x:.3f}' for x in self.para[i]])} +/- {', '.join([f'{x:.3f}' for x in self.perr[i]])} | $\chi^2$ {self.chi_square[i]:.3f}" ,
                  transform=plt.gca().transAxes, fontsize=12, fontfamily='monospace',
                  verticalalignment='top', horizontalalignment='right', color=self.default_colors[i])

    # Replace plt.title() with plt.text() to position the title inside the plot
    self.ax.text(0.02, 0.05, f'Exp Data : {self.dataName}', transform=plt.gca().transAxes,
        fontsize=12, verticalalignment='bottom', horizontalalignment='left', color='black')


    if self.legend_checkbox.isChecked() :
      fitXsecID = list(dict.fromkeys(fitXsecID))
      fitHeaders = list(dict.fromkeys(fitHeaders))

      size = len(fitHeaders) - 1
      for i , header in enumerate(fitHeaders):
        self.ax.text(0.02, 0.10 + 0.05 * size - 0.05*i, f'Fit-{fitXsecID[i]} : {header}', transform=plt.gca().transAxes,
          fontsize=12, verticalalignment='bottom', horizontalalignment='left', color='grey')


    self.ax.errorbar(self.x_exp, self.y_exp, xerr=self.x_err, yerr=self.y_err, 
                    fmt='x', label='Experimental Data', color='black', markersize = 15, elinewidth=2)


    # Plot decorator
    # Apply log scale for y-axis if selected
    if self.log_scale_checkbox.isChecked():
      self.ax.set_yscale('log')
    else:
      self.ax.set_yscale('linear')


    if self.gridline_checkbox.isChecked():
      self.ax.grid(True, which='both', linestyle='--', linewidth=0.5, color='gray')
    else:
      self.ax.grid(False)


    self.ax.set_xlabel(r'$\theta_{cm}$ [deg]')
    self.ax.set_ylabel(r'd$\sigma$/d$\Omega$ [a.u.]')
    # self.ax.legend(loc='upper right', frameon=True)

    self.ax.autoscale(enable=True, axis='x', tight=True)
    self.figure.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

    self.canvas.draw_idle()
