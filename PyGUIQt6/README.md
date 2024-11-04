
# Force install as system packages

Do it as your own risk. If you have other pyQt project, I highly NOT recommand.

 >python3 -m pip install --break-system-package PyQt6 PyQt6-WebEngine Plotly

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