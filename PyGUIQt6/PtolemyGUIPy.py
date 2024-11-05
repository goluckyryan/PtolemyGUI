#!/usr/bin/python3

import os
import platform
import subprocess
import sys
from PyQt6.QtWidgets import (
  QApplication, QMainWindow, QGridLayout, QPushButton, 
  QComboBox, QWidget, QLabel, QLineEdit, QTextEdit, QCheckBox,
  QFileDialog, QGroupBox, QVBoxLayout, QSpinBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont

from ExtractXsecPy import extract_xsec
from PlotWindow import PlotWindow
from ExWindow import ExWindow
from MatPlotLibWindow import MatPlotLibWindow

################################################## MainWindow
class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.setGeometry(100, 100, 1000, 700)
    self.setMinimumSize(400, 600)

    self.DWBAFileName = "DWBA"
    self.bashResult = ""
    self.plot_window = MatPlotLibWindow() 
    self.Ex_window = ExWindow()

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

    # Ex Group
    self.gbEx = QGroupBox("Ex")
    Ex_layout = QGridLayout()
    Ex_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbEx.setLayout(Ex_layout)

    lbName = QLabel("Isotop :")
    lbName.setAlignment(Qt.AlignmentFlag.AlignRight)

    self.leName = QLineEdit()
    self.leName.setText("12C")

    lbMaxEx = QLabel("Max Ex [MeV]:")
    lbMaxEx.setAlignment(Qt.AlignmentFlag.AlignRight)

    self.sbMaXEx = QDoubleSpinBox()
    self.sbMaXEx.setMinimum(0)
    self.sbMaXEx.setMaximum(20)
    self.sbMaXEx.setDecimals(1)
    self.sbMaXEx.setValue(10)

    buEx = QPushButton("Get & Plot Ex")
    buEx.setFixedHeight(40)
    buEx.clicked.connect(self.open_Ex_window)

    Ex_layout.addWidget(lbName, 0, 0)
    Ex_layout.addWidget(self.leName, 0, 1)
    Ex_layout.addWidget(lbMaxEx, 1, 0)
    Ex_layout.addWidget(self.sbMaXEx, 1, 1)
    Ex_layout.addWidget(buEx, 2, 0, 1, 2)

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
    font = QFont("Courier New", 10)  # You can adjust the size as needed
    self.text_edit.setFont(font)

    self.leStatus = QLineEdit("")
    self.leStatus.setReadOnly(True)

    self.LoadFileToTextBox(self.DWBAFileName)

    # Set up the layout
    layout = QGridLayout()
    # layout.addWidget(self.gbDWBA, 0, 0, 7, 1)
    layout.addWidget(self.gbDWBA, 0, 0, 5, 1)
    layout.addWidget(self.gbEx, 5, 0, 2, 1)

    layout.addWidget(self.bnOpenDWBASource, 0, 1)
    layout.addWidget(self.leFileName, 0, 2, 1, 3)
    layout.addWidget(self.bnSaveFile, 0, 5)
    layout.addWidget(self.text_edit, 1, 1, 5, 5)
    layout.addWidget(self.leStatus, 6, 1, 1, 5)

    layout.setColumnStretch(0, 1)
    layout.setColumnStretch(1, 3)

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
    
    self.SaveFile()
    
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

    isRunOK = True
    if self.chkRunPtolemy.isChecked() :
      os_name = platform.system()

      if os_name == "Linux" :
        self.BashCommand("../Cleopatra/ptolemy <" + self.DWBAFileName + ".in>" + " " + self.DWBAFileName + ".out")
      
      if os_name == "Darwin":
        self.BashCommand("../Cleopatra/ptolemy_mac <" + self.DWBAFileName + ".in>" + " " + self.DWBAFileName + ".out")

      if self.bashResult.returncode != 0 :
        isRunOK = False
        self.leStatus.setText("Ptolemy Run Error. Should check the out File.")

    if isRunOK and self.chkExtracrXsec.isChecked() and self.file_exists(self.DWBAFileName + ".out") :
      extract_xsec(self.DWBAFileName + ".out", self.cbXsec.currentIndex())
      # option = str(self.cbXsec.currentIndex())
      # self.BashCommand("../Cleopatra/ExtractXSec " + self.DWBAFileName + ".out " +  option)

    if self.chkPlot.isChecked() and self.file_exists(self.DWBAFileName + ".Xsec.txt") :
      self.open_plot_window()

  # def open_plot_window(self):
  #   if self.plot_window is None :
  #     self.plot_window = PlotWindow(self.DWBAFileName + ".Xsec.txt") 
  #     self.plot_window.show()
  #     # self.plot_window.setAttribute(Qt.WA_DeleteOnClose)  # Optional: Automatically delete when closed
  #   else:
  #     self.plot_window.read_data(self.DWBAFileName + ".Xsec.txt") 
  #     self.plot_window.plot_plotly_graph()
  #     self.plot_window.show()

  def open_plot_window(self):
    self.plot_window.read_data(self.DWBAFileName + ".Xsec.txt") 
    self.plot_window.plot_matplotlib_graph()
    self.plot_window.show()

  def open_Ex_window(self):
    self.Ex_window.GetEx(self.leName.text(), self.sbMaXEx.value())
    self.Ex_window.plot_Ex_graph()
    self.Ex_window.show()


  def closeEvent(self, event):
    if self.plot_window:
      self.plot_window.close()  # Close the PlotWindow when MainWindow closes
    if self.Ex_window:
      self.Ex_window.close()  # Close the PlotWindow when MainWindow closes
      self.Ex_window.__del__()

    event.accept()  # Accept the event to proceed with closing the main window


################################################## Main
if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())