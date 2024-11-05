
# Force install as system packages

Do it as your own risk. If you have other pyQt project, I highly NOT recommand.

 >python3 -m pip install --break-system-package PyQt6 PyQt6-WebEngine Plotly pandas matplotlib

# Use Python Virtual Enviroment

To enable virtual enviroment
 >sudo apt install python3.12-venv

create an virtual enivriment with name ptolemyGUI
 >python3 -m venv ptolemyGUI

To activate the virtual enviroment
 >source ptolemyGUI/bin/activate

To install the packages
 >pip install -r requirements.txt

To deactivate the virtual enviroment
 >deactivate

To remove the virtual enviroment, simply remove the folder
 >rm -rf ptolemyGUI

# program structure

the python script use the __Cleopatra/inFileCreator__ to generate the infile for ptolmey, and run the ptolmey, and use the ExtractXsecPy.py to output Xsec.txt. After that, the read_data() in __MatPlotLinWindow__ will read the Xsec.txt file and plot using MatPlotLib. 

There is a side program to get isotopes excitation energies from IAEA and plot the energy levels using Plotly.
