#!/usr/bin/python3

import numpy as np
from scipy.optimize import curve_fit

from ExtractXsecPy import read_DWBA
from PlotWindow import FitPlotWindow

#========================================================
class Fitting():
  def __init__(self):

    self.dataName_list = []
    self.fitOption = []
    self.expData = []

    self.dataX = []
    self.data = [] # is a 2D array
    self.headers = []

    # fit parameters for a single data set
    self.para = []
    self.para_err = []
    self.chi_squared = []

    self.plot = []

  def read_data(self,file_path):
    self.headers, self.dataX, self.data = read_DWBA(file_path)
    self.headers = self.headers[1:]

  def read_expData(self, fileName):
    self.dataName_list = []
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
            
        # Check for dataSet lines
        if line.startswith("#="):
          # If there's an existing data block, save it
          if current_data:
            self.expData.append(np.array(current_data, dtype=float))
            current_data = []
              
          # Extract dataSet Name
          dataName = line.split()[1:]
          self.dataName_list.append(" ".join(dataName))
        
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
    self.dataName_list = np.array(self.dataName_list)
    self.expData = [np.array(data) for data in self.expData]

    # Output the result
    print("=========== Number of data set:", len(self.dataName_list))
    for i in range(0, len(self.dataName_list)):
      print("-------------------------")
      print("  data Name:", self.dataName_list[i])
      print("Fit Options:", self.fitOption[i])
      print("  Data List:\n", self.expData[i])

  def FitSingleData(self, expDataID):
    print("============================================")

    # Get the number of fit options and cross-sections
    nFit = len(self.fitOption[expDataID])
    nXsec = len(self.data)

    # Extract experimental data (x, y, and errors)
    x_exp = self.expData[expDataID][:, 0]  # x positions
    y_exp = self.expData[expDataID][:, 2]  # y positions
    y_err = self.expData[expDataID][:, 3]  # y errors

    self.para = []
    self.para_err = []
    self.chi_squared = []

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

      self.para.append(popt)
      perr = np.sqrt(np.diag(pcov))# Standard deviation of the parameters
      self.para_err.append(perr)  

      # Get the fitted model values
      y_fit = fit_func(x_exp, *popt)
      residuals = y_exp - y_fit
      self.chi_squared.append(np.sum((residuals / y_err) ** 2))

      print(f"Fitted scale for fit {k}: {', '.join([f'{x:.3f}' for x in popt])} +/- {', '.join([f'{x:.3f}' for x in perr])} | Chi^2 : {self.chi_squared[-1]:.4f}")
      # print(f"Fitted scale for fit {k}: {popt} +/- {perr} | Chi^2 : {chi_squared[-1]:.4f}")

    return self.para, self.para_err, self.chi_squared


  def plot_fits(self):

    self.plot = []

    for k , dN in enumerate(self.dataName_list):
      self.FitSingleData(k)
      self.plot.append( FitPlotWindow(f"Data-{k}"))
      self.plot[-1].set_data(k, self.expData, self.fitOption, dN, 
                             self.dataX, self.data, self.headers,
                             self.para, self.para_err, self.chi_squared)
      self.plot[-1].plot_Fit()
      self.plot[-1].show()

  def close_plots(self):
    for plot in self.plot:
      plot.close()
