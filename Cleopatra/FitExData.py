#!/usr/bin/python3

import numpy as np
from scipy.optimize import curve_fit

from PyQt6.QtWidgets import (
  QVBoxLayout, QWidget
)

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

from ExtractXsecPy import read_DWBA

default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

class FitPlotWidget(QWidget):
  def __init__(self, figure):
    super().__init__()

    self.setWindowTitle("Fit Plot")
    self.resize(800, 600)

    self.canvas = FigureCanvas(figure)
    self.toolbar = NavigationToolbar(self.canvas, self)

    layout = QVBoxLayout(self)
    layout.addWidget(self.toolbar)
    layout.addWidget(self.canvas)

    self.setLayout(layout)


class Fitting():
  def __init__(self):

    self.ExList = []
    self.fitOption = []
    self.expData = []

    self.dataX = []
    self.data = [] # is a 2D array
    self.headers = []

  def read_data(self,file_path):
    self.headers, self.dataX, self.data = read_DWBA(file_path)

    print(self.headers)

  def read_expData(self, fileName):
    self.ExList = []
    self.fitOption = []
    self.expData = []

    current_data = []

    with open(fileName, "r") as file:
      for line in file:
        line = line.strip()
        
        if not line:
          continue

        if line.startswith("$"):
          continue
            
        # Check for excitation energy lines
        if line.startswith("#="):
          # If there's an existing data block, save it
          if current_data:
            self.expData.append(np.array(current_data, dtype=float))
            current_data = []
              
          # Extract excitation energy
          Ex_energy = line.split()[1]
          self.ExList.append(float(Ex_energy))
        
        # Check for fit option lines
        elif line.startswith("fit"):
          # Add fit option parameters
          fit_params = [x.strip(',') for x in line.split()[1:]]
          self.fitOption.append(fit_params)
        
        # Parse data lines
        elif not line.startswith("#"):
          values = [float(x) for x in line.split()]
          current_data.append(values)

      # Append the last block
      if current_data:
        self.expData.append(np.array(current_data, dtype=float))

    # Convert to numpy arrays
    self.ExList = np.array(self.ExList)
    self.expData = [np.array(data) for data in self.expData]

    # Output the result
    print("=========== Number of data set:", len(self.ExList))
    for i in range(0, len(self.ExList)):
      print("-------------------------")
      print("     ExList:", self.ExList[i])
      print("Fit Options:", self.fitOption[i])
      print("  Data List:\n", self.expData[i])


  def FitData(self):

    figure_list = []

    for expDataID in range(len(self.expData)):
      
      print("============================================")

      # Get the number of fit options and cross-sections
      nFit = len(self.fitOption[expDataID])
      nXsec = len(self.data)

      # Extract experimental data (x, y, and errors)
      x_exp = self.expData[expDataID][:, 0]  # x positions
      x_err = self.expData[expDataID][:, 1]  # x uncertainties (errors)
      y_exp = self.expData[expDataID][:, 2]  # y positions
      y_err = self.expData[expDataID][:, 3]  # y errors

      fitTheory = []
      fitTheory_lower = []
      fitTheory_upper = []

      para = []
      perr = []
      chi_squared = []

      for k in range(nFit):
          # Get the cross-section IDs for the current fit option and strip extra spaces
          xsecIDStr = self.fitOption[expDataID][k].strip()
          xsecID = [int(part.strip()) for part in xsecIDStr.split('+')] if '+' in xsecIDStr else [int(xsecIDStr)]

          # Ensure all cross-section IDs are valid
          processFlag = True
          for id in range(len(xsecID)):
              if xsecID[id] >= nXsec:
                  print(f"Error: Requested Xsec-{xsecID[id]} exceeds the number of available cross-sections ({nXsec})")
                  processFlag = False
          
          if processFlag == False :
            continue

          # Define the fitting function: a weighted sum of the selected data
          def fit_func(x, *scale):
              y = np.zeros_like(x)
              for p, id in enumerate(xsecID):
                  y += scale[p] * np.interp(x, self.dataX, self.data[id])
              return y


          lower_bounds = [1e-6] * len(xsecID)  # Setting a small positive lower bound
          upper_bounds = [np.inf] * len(xsecID)  # No upper bound

          # Perform curve fitting using the fit_func and experimental data with y-errors as weights
          popt, pcov = curve_fit(fit_func, x_exp, y_exp, sigma=y_err, absolute_sigma=True, 
                                p0=np.ones(len(xsecID)), # Initial guess for scale parameters
                                bounds=(lower_bounds, upper_bounds)) 

          para.append(popt)
          perr.append(np.sqrt(np.diag(pcov)))  # Standard deviation of the parameters

          # Get the fitted model values
          y_fit = fit_func(x_exp, *popt)
          residuals = y_exp - y_fit
          chi_squared.append(np.sum((residuals / y_err) ** 2))

          print(f"Fitted scale for fit {k}: {', '.join([f'{x:.3f}' for x in popt])} +/- {', '.join([f'{x:.3f}' for x in perr[-1]])} | Chi^2 : {chi_squared[-1]:.4f}")
          # print(f"Fitted scale for fit {k}: {popt} +/- {perr} | Chi^2 : {chi_squared[-1]:.4f}")

          # Append the theoretical fit for this fit option
          fitTheory.append(np.zeros_like(self.dataX))
          for p, id in enumerate(xsecID):
              fitTheory[k] += popt[p] * np.interp(self.dataX, self.dataX, self.data[id])

          # Optionally, you can plot the uncertainty as shaded regions (confidence intervals)
          # Create the upper and lower bounds of the theoretical model with uncertainties
          fitTheory_upper.append(np.zeros_like(self.dataX))
          fitTheory_lower.append(np.zeros_like(self.dataX))
          
          for p, id in enumerate(xsecID):
              fitTheory_upper[k] += (popt[p] + perr[p]) * np.interp(self.dataX, self.dataX, self.data[id])
              fitTheory_lower[k] += (popt[p] - perr[p]) * np.interp(self.dataX, self.dataX, self.data[id])

      fig = plt.figure()
      figure_list.append(fig)

      # Plot results
      plt.errorbar(x_exp, y_exp, xerr=x_err, yerr=y_err, fmt='x', label='Experimental Data', color='black', markersize = 15, elinewidth=2)

      # Plot all fit theories
      for i, fit in enumerate(fitTheory):
        plt.plot(self.dataX, fit, label=f'Chi2:{chi_squared[i]:.3f} | Xsec:{self.fitOption[expDataID][i]}')
        plt.fill_between(self.dataX, fitTheory_lower[i], fitTheory_upper[i], alpha=0.2)

      # Customize plot
      plt.xlabel('Angle_CM [deg]')
      plt.ylabel('X-Sec [a.u.]')
      plt.legend()
      plt.autoscale(enable=True, axis='x', tight=True)
      plt.tight_layout()
      plt.yscale('log')

      # Replace plt.title() with plt.text() to position the title inside the plot
      plt.text(0.05, 0.05, f'Fit for Exp Data : {self.ExList[expDataID]} MeV', transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='bottom', horizontalalignment='left', color='black')

      for i, _ in enumerate(para):
        plt.text(0.05, 0.1 + 0.05*i, f"Xsec-{self.fitOption[expDataID][i].strip()}: {', '.join([f'{x:.3f}' for x in para[i]])} +/- {', '.join([f'{x:.3f}' for x in perr[i]])}" , transform=plt.gca().transAxes,
           fontsize=12, verticalalignment='bottom', horizontalalignment='left', color=default_colors[i])


      
    return figure_list