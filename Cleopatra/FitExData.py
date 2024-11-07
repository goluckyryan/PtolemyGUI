#!/usr/bin/python3

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from ExtractXsecPy import read_DWBA

class Fitting:
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
            
        # Check for excitation energy lines
        if line.startswith("#======================"):
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

    # Set initial offset values
    x_offset, y_offset = 2000, 100

    figure = []

    for expDataID in range(len(self.expData)):
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

          perr = np.sqrt(np.diag(pcov))  # Standard deviation of the parameters
          print(f"Fitted scale for fit {k+1}: {popt} +/- {perr}")

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
      figure.append(fig)

      # Plot results
      plt.errorbar(x_exp, y_exp, xerr=x_err, yerr=y_err, fmt='o', label='Experimental Data', color='blue')

      # Plot all fit theories
      for i, fit in enumerate(fitTheory):
        plt.plot(self.dataX, fit, label=f'Xsec:{self.fitOption[expDataID][i]} Fit')
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

      manager = plt.get_current_fig_manager()
      # manager.window.wm_geometry(f"+{x_offset}+{y_offset}")
      manager.set_window_title(f"Exp Data : {self.ExList[expDataID]} MeV")
      
      x_offset += 100
      y_offset += 100

    plt.show()