#!/usr/bin/python3

import os
import platform
import subprocess
import sys
from functools import partial
from PyQt6.QtWidgets import (
  QApplication, QMainWindow, QGridLayout, QPushButton, 
  QComboBox, QWidget, QLabel, QLineEdit, QTextEdit, QCheckBox,
  QFileDialog, QGroupBox, QVBoxLayout, QSpinBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt, QUrl
from PyQt6.QtGui import QFont
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
import tempfile

class PlotWindow(QWidget):
  def __init__(self, XsecFile):
    super().__init__()

    self.setWindowTitle("DWBA Plot")
    self.setGeometry(100, 100, 800, 600)

    self.x = []
    self.data = []
    self.headers = []
    self.x, self.data, self.headers = self.read_data(XsecFile)

    layout = QVBoxLayout(self)
    self.web_view = QWebEngineView()
    layout.addWidget(self.web_view)

    self.plot_plotly_graph()

  def read_data(self,file_path):
    x = []  # List for the first column
    data = [] # 2D list for other columns
    headers = []  # List to store headers

    with open(file_path, 'r') as file:
      header_found = False  # Flag to indicate if the header has been found
      for line in file:
        # Skip lines that start with '#' and empty lines
        if line.startswith('#') or not line.strip():
          continue
        
        if not header_found:
          parts = line.split('ELab')
          elab_parts = [parts[0]]  # Start with the first part
          for part in parts[1:]:
              elab_parts.append('ELab' + part)  # Prepend 'ELab' to each subsequent part

          headers = elab_parts  # Use the split parts as headers
          header_found = True  # Set the flag to True to skip this block in future iterations
          print(f"ELab parts found: {elab_parts}")  # Print or process this as needed
          continue
        
        # Split the line by whitespace
        parts = line.split()
        if len(parts) > 0:  # Make sure there is at least one column
          x.append(float(parts[0]))  # First column
          # Append the rest of the columns to data
          if len(data) == 0:
            # Initialize the data array with the right number of sublists
            data = [[] for _ in range(len(parts) - 1)]
          for i in range(len(parts) - 1):
            data[i].append(float(parts[i + 1]))  # Rest of the columns


      print(headers)
    return x, data, headers

  def plot_plotly_graph(self):
    # Create a Plotly figure
    fig = go.Figure()

    # Add traces for each column in data against x
    for i, y in enumerate(self.data):
        fig.add_trace(go.Scatter(x=self.x, y=y, mode='lines+markers', name=self.headers[i + 1]))  # Use headers for names

    # Update layout for better presentation
    fig.update_layout(
        xaxis_title="Angle_CM [Deg]",
        yaxis_title="Values",
        template="plotly"
    )

    # Save the plot as an HTML file in a temporary location
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
        fig.write_html(tmp_file.name)
        html_file = tmp_file.name

    # Load the HTML file in QWebEngineView
    self.web_view.setUrl(QUrl.fromLocalFile(html_file))

class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.setGeometry(100, 100, 1000, 700)
    self.setMinimumSize(400, 600)

    self.DWBAFileName = "DWBA"
    self.bashResult = ""
    self.plot_window = None

    # Set up Group Box for DWBA Control
    self.gbDWBA = QGroupBox("DWBA")
    group_layout = QGridLayout()
    group_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbDWBA.setLayout(group_layout)

    self.bnOpenDWBA = QPushButton("Open DWBA")
    self.bnOpenDWBA.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName))
    self.bnOpenInFile = QPushButton("Open *.in File")
    self.bnOpenInFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".in"))
    self.bnOpenOutFile = QPushButton("Open *.out File")
    self.bnOpenOutFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".out"))
    self.bnOpenXsecFile = QPushButton("Open X-sec File")
    self.bnOpenXsecFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".Xsec.txt"))

    lbAngMin = QLabel("angMin")
    lbAngMax = QLabel("angMax")
    lbAngSize = QLabel("angSize")
    self.sbAngMin = QSpinBox()
    self.sbAngMin.setValue(0)
    self.sbAngMin.setMinimum(0)
    self.sbAngMin.setMaximum(180)
    self.sbAngMax = QSpinBox()
    self.sbAngMax.setValue(60)
    self.sbAngMax.setMinimum(0)
    self.sbAngMax.setMaximum(180)
    self.sbAngSize = QDoubleSpinBox()
    self.sbAngSize.setValue(1)
    self.sbAngSize.setMinimum(0.1)
    self.sbAngSize.setMaximum(10)
    self.sbAngSize.setSingleStep(0.5)

    self.chkCreateInFile = QCheckBox("Create InFile")
    self.chkCreateInFile.setChecked(True)
    self.chkRunPtolemy = QCheckBox("Run Ptolemy")
    self.chkRunPtolemy.setChecked(True)
    self.chkExtracrXsec = QCheckBox("Extract Xsec")
    self.chkExtracrXsec.setChecked(True)
    self.chkExtracrXsec.stateChanged.connect(self.OnOffXsecOption)

    self.cbXsec = QComboBox()
    self.cbXsec.addItem("XSec")
    self.cbXsec.addItem("Ratio to Ruth.")
    self.cbXsec.addItem("Ruth.")

    self.chkPlot = QCheckBox("Plot")
    self.chkPlot.setChecked(True)

    self.bnCalDWBA = QPushButton("Calculate DWBA")
    self.bnCalDWBA.setFixedHeight(50)
    self.bnCalDWBA.clicked.connect(self.CalDWBA)

    group_layout.addWidget(self.bnOpenDWBA, 0, 0, 1, 2)
    group_layout.addWidget(self.bnOpenInFile, 1, 0, 1, 2)
    group_layout.addWidget(self.bnOpenOutFile, 2, 0, 1, 2)
    group_layout.addWidget(self.bnOpenXsecFile, 3, 0, 1, 2)

    group_layout.addWidget(lbAngMin, 4, 0)
    group_layout.addWidget(self.sbAngMin, 4, 1)
    group_layout.addWidget(lbAngMax, 5, 0)
    group_layout.addWidget(self.sbAngMax, 5, 1)
    group_layout.addWidget(lbAngSize, 6, 0)
    group_layout.addWidget(self.sbAngSize, 6, 1)

    group_layout.addWidget(self.chkCreateInFile, 7, 0, 1, 2)
    group_layout.addWidget(self.chkRunPtolemy, 8, 0, 1, 2)
    group_layout.addWidget(self.chkExtracrXsec, 9, 0, 1, 2)

    group_layout.addWidget(self.cbXsec, 10, 0, 1, 2)
    group_layout.addWidget(self.chkPlot, 11, 0, 1, 2)

    group_layout.addWidget(self.bnCalDWBA, 12, 0, 1, 2)

    # Set up the Right Side

    self.bnOpenDWBASource = QPushButton("Open DWBA Source")
    self.bnOpenDWBASource.clicked.connect(self.OpenDWBASourceFile)

    self.leFileName = QLineEdit("")
    self.leFileName.setReadOnly(True)
    self.leFileName.setText(self.DWBAFileName)

    self.bnSaveFile = QPushButton("Save File")
    self.bnSaveFile.clicked.connect(self.SaveFile)

    self.text_edit = QTextEdit()
    self.text_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
    font = QFont("Courier New", 8)  # You can adjust the size as needed
    self.text_edit.setFont(font)

    self.leStatus = QLineEdit("")
    self.leStatus.setReadOnly(True)

    self.LoadFileToTextBox(self.DWBAFileName)

    # Set up the layout
    layout = QGridLayout()
    layout.addWidget(self.gbDWBA, 0, 0, 7, 1)

    layout.addWidget(self.bnOpenDWBASource, 0, 1)
    layout.addWidget(self.leFileName, 0, 2, 1, 3)
    layout.addWidget(self.bnSaveFile, 0, 5)
    layout.addWidget(self.text_edit, 1, 1, 5, 5)
    layout.addWidget(self.leStatus, 6, 1, 1, 5)

    # Set up the container and layout
    container = QWidget()
    container.setLayout(layout)
    self.setCentralWidget(container)

  ####################################### methods
  def OnOffXsecOption(self):
    self.cbXsec.setEnabled(self.chkExtracrXsec.isChecked())

  def OpenDWBASourceFile(self):
    file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*)")        
    if file_path:
      self.DWBAFileName = file_path
      self.leFileName.setText(self.DWBAFileName)
      self.LoadFileToTextBox(self.DWBAFileName)

  def LoadFileToTextBox(self, fileName):
    # print(fileName)
    try:
      with open(fileName, 'r') as file:
        content = file.read()
        self.text_edit.setText(content)
        self.leStatus.setText(f"Loaded file : {fileName}")
        self.leFileName.setText(fileName)
    except Exception as e:
      self.text_edit.setText(f"Failed to load file:\n{e}")
      self.leStatus.setText(f"Failed to load file:\n{e}")
  
  def SaveFile(self):
    file_path = self.leFileName.text()
    with open(file_path, 'w') as file:
      file.write(self.text_edit.toPlainText())
      self.leStatus.setText(f"File saved to: {file_path}")

  def BashCommand(self, cmd):
    print("Bash Command : |" + cmd + "|")
    self.bashResult = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    print("Output:", self.bashResult.stdout)
    print("Error:", self.bashResult.stderr)
    print("Return Code:", self.bashResult.returncode)

  def file_exists(self,file_path):
    return os.path.exists(file_path) and os.path.isfile(file_path)
  
  def CalDWBA(self):
    
    self.BashCommand("cd ../Cleopatra; make;cd ../PyGUIQt6")

    print(" Is Create InFile : " + str(self.chkCreateInFile.isChecked() ))
    print("   Is Run Ptolemy : " + str(self.chkRunPtolemy.isChecked() ))
    print("  Is Extract XSec : " + str(self.chkExtracrXsec.isChecked() ))
    print("          Is Plot : " + str(self.chkPlot.isChecked() ))

    if self.chkCreateInFile.isChecked() :
      aMin  = " " + str(self.sbAngMin.value())
      aMax  = " " + str(self.sbAngMax.value())
      aSize = " " + str(self.sbAngSize.value())

      self.BashCommand("../Cleopatra/InFileCreator " +  self.DWBAFileName + aMin + aMax + aSize)

    isRunOK = False
    if self.chkRunPtolemy.isChecked() :
      os_name = platform.system()

      if os_name == "Linux" :
        self.BashCommand("../Cleopatra/ptolemy <" + self.DWBAFileName + ".in>" + " " + self.DWBAFileName + ".out")
      
      if os_name == "Darwin":
        self.BashCommand("../Cleopatra/ptolemy_mac <" + self.DWBAFileName + ".in>" + " " + self.DWBAFileName + ".out")

      if self.bashResult.returncode == 0 :
        isRunOK = True
      else:
        self.leStatus.setText("Ptolemy Run Error. Should check the out File.")

    if isRunOK and self.chkExtracrXsec.isChecked() and self.file_exists(self.DWBAFileName + ".out") :
      option = str(self.cbXsec.currentIndex())
      self.BashCommand("../Cleopatra/ExtractXSec " + self.DWBAFileName + ".out " +  option)

    ### Plot ##
    if self.chkPlot.isChecked() and self.file_exists(self.DWBAFileName + ".Xsec.txt") :
      print("dasdsadsa")
      self.open_plot_window()

  def open_plot_window(self):
    print("dsadasdsadafkal;k;lkl;kdasdsadsa")
    if self.plot_window is None :
      self.plot_window = PlotWindow(self.DWBAFileName + ".Xsec.txt") 
      self.plot_window.show()
      self.plot_window.setAttribute(Qt.WA_DeleteOnClose)  # Optional: Automatically delete when closed
    else:
      self.plot_window.read_data(self.DWBAFileName + ".Xsec.txt") 
      self.plot_window.plot_plotly_graph()

  def closeEvent(self, event):
    if self.plot_window:
      self.plot_window.close()  # Close the PlotWindow when MainWindow closes
    event.accept()  # Accept the event to proceed with closing the main window


################################################## Main
if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())