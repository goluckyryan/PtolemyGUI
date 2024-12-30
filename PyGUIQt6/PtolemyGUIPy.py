#!/usr/bin/env python3

import os
import platform
import subprocess
import sys
from PyQt6.QtWidgets import (
  QApplication, QMainWindow, QGridLayout, QPushButton, 
  QComboBox, QWidget, QLabel, QLineEdit, QCheckBox,
  QFileDialog, QGroupBox, QVBoxLayout, QSpinBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))

from CustomTextEdit import CustomTextEdit
from ExtractXsecPy import extract_xsec
from ExWindow import ExWindow
from MatPlotLibWindow import MatPlotLibWindow
from FitExData import Fitting

################################################## MainWindow
class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.resize(1000, 800)
    self.setMinimumSize(1000, 800)

    self.lastDWBARecord = "lastDWBA.txt"
    self.DWBAFileName = ""
    self.ExpDataFileName = ""
    self.LoadLastOpenDWBASource()
    self.bashResult = ""
    self.plot_window = MatPlotLibWindow() 
    self.Ex_window = ExWindow()
    self.fitting = Fitting()
    self.fitCanvas = []

    # Set up Group Box for DWBA Control
    self.gbDWBA = QGroupBox("DWBA")
    group_layout = QGridLayout()
    group_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbDWBA.setLayout(group_layout)

    self.bnOpenDWBA = QPushButton("Open DWBA")
    self.bnOpenDWBA.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName, True))
    self.bnOpenInFile = QPushButton("Open *.in File")
    self.bnOpenInFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".in"))
    self.bnOpenOutFile = QPushButton("Open *.out File")
    self.bnOpenOutFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".out"))
    self.bnOpenXsecFile = QPushButton("Open X-sec File")
    self.bnOpenXsecFile.clicked.connect(lambda: self.LoadFileToTextBox(self.DWBAFileName + ".Xsec.txt"))

    self.bnDeleteFiles = QPushButton("Delete in/out/Xsec files")
    self.bnDeleteFiles.clicked.connect(self.DeleteinOutXsecFiles)

    lbAngMin = QLabel("angMin :")
    lbAngMin.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)
    lbAngMax = QLabel("angMax :")
    lbAngMax.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)
    lbAngSize = QLabel("angSize :")
    lbAngSize.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)
    self.sbAngMin = QSpinBox()
    self.sbAngMin.setValue(0)
    self.sbAngMin.setMinimum(0)
    self.sbAngMin.setMaximum(180)
    self.sbAngMax = QSpinBox()
    self.sbAngMax.setValue(60)
    self.sbAngMax.setMinimum(0)
    self.sbAngMax.setMaximum(180)
    self.sbAngSize = QDoubleSpinBox()
    self.sbAngSize.setValue(0.2)
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
    group_layout.addWidget(self.bnDeleteFiles, 4, 0, 1, 2)

    group_layout.addWidget(lbAngMin, 5, 0)
    group_layout.addWidget(self.sbAngMin, 5, 1)
    group_layout.addWidget(lbAngMax, 6, 0)
    group_layout.addWidget(self.sbAngMax, 6, 1)
    group_layout.addWidget(lbAngSize, 7, 0)
    group_layout.addWidget(self.sbAngSize, 7, 1)

    group_layout.addWidget(self.chkCreateInFile, 8, 0, 1, 2)
    group_layout.addWidget(self.chkRunPtolemy, 9, 0, 1, 2)
    group_layout.addWidget(self.chkExtracrXsec, 10, 0, 1, 2)

    group_layout.addWidget(self.cbXsec, 11, 0, 1, 2)
    group_layout.addWidget(self.chkPlot, 12, 0, 1, 2)

    group_layout.addWidget(self.bnCalDWBA, 13, 0, 1, 2)

    # Ex Group
    self.gbEx = QGroupBox("Ex")
    Ex_layout = QGridLayout()
    Ex_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbEx.setLayout(Ex_layout)

    lbName = QLabel("Isotop :")
    lbName.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)

    self.leName = QLineEdit()
    self.leName.setText("12C")

    lbMaxEx = QLabel("Max Ex [MeV]:")
    lbMaxEx.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignCenter)

    self.sbMaXEx = QDoubleSpinBox()
    self.sbMaXEx.setMinimum(0)
    self.sbMaXEx.setMaximum(20)
    self.sbMaXEx.setDecimals(1)
    self.sbMaXEx.setValue(10)

    buEx = QPushButton("Get And Plot Ex")
    buEx.setFixedHeight(40)
    buEx.clicked.connect(self.open_Ex_window)

    Ex_layout.addWidget(lbName, 0, 0)
    Ex_layout.addWidget(self.leName, 0, 1)
    Ex_layout.addWidget(lbMaxEx, 1, 0)
    Ex_layout.addWidget(self.sbMaXEx, 1, 1)
    Ex_layout.addWidget(buEx, 2, 0, 1, 2)

    # ExpData Group
    self.gbExpFit = QGroupBox("Exp Data Fit")
    fit_layout = QVBoxLayout()
    fit_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbExpFit.setLayout(fit_layout)

    bnOpenExpData = QPushButton("Open ExpData")
    bnOpenExpData.clicked.connect(self.LoadExpDataToTextBox)

    bnFit = QPushButton("Fit")
    bnFit.clicked.connect(self.fitData)

    fit_layout.addWidget(bnOpenExpData)
    fit_layout.addWidget(bnFit)

    # Set up the Right Side

    self.bnOpenDWBASource = QPushButton("Open DWBA Source")
    self.bnOpenDWBASource.clicked.connect(self.OpenDWBASourceFile)

    self.leFileName = QLineEdit("")
    self.leFileName.setReadOnly(True)
    self.leFileName.setText(self.DWBAFileName)

    self.bnSaveFile = QPushButton("Save File")
    self.bnSaveFile.clicked.connect(self.SaveFile)

    self.bnOpenExpData = QPushButton("Open Exp Data")
    self.bnOpenExpData.clicked.connect(self.OpenExpDataFile)

    self.leExpDataFileName = QLineEdit("")
    self.leExpDataFileName.setReadOnly(True)
    self.leExpDataFileName.setText(self.ExpDataFileName)

    self.bnSaveExpDataFile = QPushButton("Save Exp Data File")
    self.bnSaveExpDataFile.clicked.connect(self.SaveExpDataFile)

    self.text_edit = CustomTextEdit(self)

    self.leStatus = QLineEdit("")
    self.leStatus.setReadOnly(True)

    self.LoadFileToTextBox(self.DWBAFileName, True)

    self.bnSaveFile.setEnabled(True)
    self.bnSaveExpDataFile.setEnabled(False)

    # Set up the layout
    layout = QGridLayout()
    # layout.addWidget(self.gbDWBA, 0, 0, 7, 1)
    layout.addWidget(self.gbDWBA, 0, 0, 5, 1)
    layout.addWidget(self.gbEx, 5, 0, 2, 1)
    layout.addWidget(self.gbExpFit, 7, 0, 2, 1)

    layout.addWidget(self.bnOpenDWBASource, 0, 1)
    layout.addWidget(self.leFileName, 0, 2, 1, 3)
    layout.addWidget(self.bnSaveFile, 0, 5)

    layout.addWidget(self.bnOpenExpData, 1, 1)
    layout.addWidget(self.leExpDataFileName, 1, 2, 1, 3)
    layout.addWidget(self.bnSaveExpDataFile, 1, 5)

    layout.addWidget(self.text_edit, 2, 1, 6, 5)
    layout.addWidget(self.leStatus, 8, 1, 1, 5)

    for i in range(layout.columnCount()) :
      layout.setColumnStretch(i, 1)
  
    # Set up the container and layout
    container = QWidget()
    container.setLayout(layout)
    self.setCentralWidget(container)

    self.text_edit.setFocus()

  ####################################### methods
  def LoadLastOpenDWBASource(self):
    try :
      with open(self.lastDWBARecord, 'r') as file:
        self.DWBAFileName = file.readline().strip()
        self.ExpDataFileName = file.readline().strip()
    except:
      self.DWBAFileName = "DWBA"
      self.ExpDataFileName = ""

  def SaveLastOpenDWBASource(self):
    with open(self.lastDWBARecord, 'w') as file:
      file.write(self.DWBAFileName)
      file.write("\n")
      file.write(self.ExpDataFileName)

  def OnOffXsecOption(self):
    self.cbXsec.setEnabled(self.chkExtracrXsec.isChecked())

  def OpenDWBASourceFile(self):
    file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*)")        
    if file_path:
      self.DWBAFileName = file_path
      self.leFileName.setText(self.DWBAFileName)
      self.LoadFileToTextBox(self.DWBAFileName)
      self.SaveLastOpenDWBASource()
      self.bnSaveExpDataFile.setEnabled(False)
      self.bnSaveFile.setEnabled(True)

  def OpenExpDataFile(self):
    file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "Text File (*.txt)")        
    if file_path:
      self.ExpDataFileName = file_path
      self.leExpDataFileName.setText(self.ExpDataFileName)
      self.LoadFileToTextBox(self.ExpDataFileName)
      self.bnSaveExpDataFile.setEnabled(True)
      self.bnSaveFile.setEnabled(False)
      self.SaveLastOpenDWBASource()

  def LoadExpDataToTextBox(self):
    if self.ExpDataFileName == "" :
      self.text_edit.clear()
      self.text_edit.append("$<-- for comment line")
      self.text_edit.append("$No expData found, this is a template")
      self.text_edit.append("$line start with '#=' starts a data set")
      self.text_edit.append("#============== state")
      self.text_edit.append("$line started with 'fit' indicate which DWBA Xsec to be fitted. ")
      self.text_edit.append("$0 for the fist one, 0+1 fit both 0 and 1.")
      self.text_edit.append("fit 0, 0+1")
      self.text_edit.append("$angle_deg  ang_err   count   count_err")
      self.text_edit.append("10          1         100     10")
      self.text_edit.append("20          1         200     14")
      self.text_edit.append("30          1         80      7")
    else:
      self.LoadFileToTextBox(self.ExpDataFileName)
      self.leFileName.setText(self.DWBAFileName)

    self.bnSaveExpDataFile.setEnabled(True)
    self.bnSaveFile.setEnabled(False)

  def LoadFileToTextBox(self, fileName, moveToButton = False):    
    self.bnSaveExpDataFile.setEnabled(False)
    self.bnSaveFile.setEnabled(True)

    # print(fileName)
    try:
      with open(fileName, 'r') as file:
        content = file.read()
        self.text_edit.setText(content)
        if moveToButton :
          cursor = self.text_edit.textCursor()
          cursor.movePosition(cursor.MoveOperation.End)
          cursor.movePosition(cursor.MoveOperation.StartOfBlock)
          self.text_edit.setTextCursor(cursor)

        self.leStatus.setText(f"Loaded file : {fileName}")
        self.leFileName.setText(fileName)
    except Exception as e:
      self.text_edit.setText(f"Failed to load file:\n{e}")
      self.leStatus.setText(f"Failed to load file:\n{e}")
    
  def SaveFile(self):
    if self.bnSaveFile.isEnabled() :
      file_path = self.leFileName.text()
      with open(file_path, 'w') as file:
        file.write(self.text_edit.toPlainText())
        self.leStatus.setText(f"File saved to: {file_path}")
  
  def SaveExpDataFile(self):
    if self.bnSaveExpDataFile.isEnabled() :
      file_path = self.leExpDataFileName.text()
      if file_path == "" :
        file_path, _ = QFileDialog.getSaveFileName(self, "Save File", "", "Text Files (*.txt)")
        if not file_path.lower().endswith(".txt"):
          file_path += ".txt"  
        self.leExpDataFileName.setText(file_path)
        self.ExpDataFileName = file_path

      with open(file_path, 'w') as file:
        file.write(self.text_edit.toPlainText())
        self.leStatus.setText(f"File saved to: {file_path}")
      
      self.SaveLastOpenDWBASource()

  def DeleteinOutXsecFiles(self):
    if os.path.exists(self.DWBAFileName + ".in"):
      os.remove(self.DWBAFileName + ".in")
    if os.path.exists(self.DWBAFileName + ".out"):
      os.remove(self.DWBAFileName + ".out")
    if os.path.exists(self.DWBAFileName + ".Xsec.txt"):
      os.remove(self.DWBAFileName + ".Xsec.txt")
    self.leStatus.setText("Deleted " + self.DWBAFileName + ".in/.out/.Xsec.txt files")

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

    if isRunOK and self.chkPlot.isChecked() and self.file_exists(self.DWBAFileName + ".Xsec.txt") :
      self.open_plot_window()

  def open_plot_window(self):
    self.plot_window.read_data(self.DWBAFileName + ".Xsec.txt") 
    self.plot_window.plot_matplotlib_graph()
    self.plot_window.show()

  def open_Ex_window(self):
    self.Ex_window.GetEx(self.leName.text(), self.sbMaXEx.value())
    if self.sbMaXEx.value() > 0 :
      self.Ex_window.plot_Ex_graph()
      self.Ex_window.show()

  def fitData(self):
    self.SaveExpDataFile()

    self.fitting.read_expData(self.ExpDataFileName)
    self.fitting.read_data(self.DWBAFileName + ".Xsec.txt")
    self.fitting.plot_fits()

  def closeEvent(self, event):
    if self.plot_window:
      self.plot_window.close()  # Close the PlotWindow when MainWindow closes
    if self.Ex_window:
      self.Ex_window.close()  # Close the PlotWindow when MainWindow closes
      self.Ex_window.__del__()
    
    self.fitting.close_plots()
    print("============== Bye Bye ========== ")
    event.accept()  # Accept the event to proceed with closing the main window


################################################## Main
if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())