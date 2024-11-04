#!/usr/bin/python3

import os
import datetime
import csv
import subprocess
import sys
import time
from functools import partial
from PyQt6.QtWidgets import (
  QApplication, QMainWindow, QGridLayout, QPushButton, 
  QComboBox, QWidget, QLabel, QLineEdit, QTextEdit, QCheckBox,
  QFileDialog, QGroupBox, QVBoxLayout, QSpinBox, QDoubleSpinBox
)
from PyQt6.QtCore import Qt, QPoint
from PyQt6.QtGui import QFont


class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.setGeometry(100, 100, 1000, 700)
    self.setMinimumSize(400, 600)

    self.DWBAFileName = "../DWBA"

    # Set up Group Box for DWBA Control
    self.gbDWBA = QGroupBox("DWBA")
    group_layout = QGridLayout()
    group_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    self.gbDWBA.setLayout(group_layout)

    self.bnOpenDWBA = QPushButton("Open DWBA Source")
    self.bnOpenDWBA.clicked.connect(self.OpenDWBASourceFile)

    self.bnOpenInFile = QPushButton("Open *.in File")
    self.bnOpenInFile.clicked.connect(partial(self.LoadFileToTextBox, self.DWBAFileName + ".in"))
    self.bnOpenOutFile = QPushButton("Open *.out File")
    self.bnOpenInFile.clicked.connect(partial(self.LoadFileToTextBox, self.DWBAFileName + ".out"))
    self.bnOpenXsecFile = QPushButton("Open X-sec File")
    self.bnOpenInFile.clicked.connect(partial(self.LoadFileToTextBox, self.DWBAFileName + ".txt"))

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
    self.chkRunPtolemy = QCheckBox("Run Ptolemy")
    self.chkExtracrXsec = QCheckBox("Extract Xsec")
    self.chkPlot = QCheckBox("Plot")

    self.cbXsec = QComboBox()
    self.cbXsec.addItem("XSec")
    self.cbXsec.addItem("Ratio to Ruth.")
    self.cbXsec.addItem("Ruth.")

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

    group_layout.addWidget(self.chkPlot, 10, 0)
    group_layout.addWidget(self.cbXsec, 10, 1)

    group_layout.addWidget(self.bnCalDWBA, 11, 0, 1, 2)

    # Set up the Right Side

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
    layout.addWidget(self.leFileName, 0, 1, 1, 4)
    layout.addWidget(self.bnSaveFile, 0, 5)
    layout.addWidget(self.text_edit, 1, 1, 5, 5)
    layout.addWidget(self.leStatus, 6, 1, 1, 5)

    # Set up the container and layout
    container = QWidget()
    container.setLayout(layout)
    self.setCentralWidget(container)

  ####################################### methods
  def OpenDWBASourceFile(self):
    file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*)")        
    if file_path:
      self.DWBAFileName = file_path
      self.leFileName.setText(self.DWBAFileName)
      self.LoadFileToTextBox(self.DWBAFileName)

  def LoadFileToTextBox(self, fileName):
    try:
      with open(fileName, 'r') as file:
        content = file.read()
        self.text_edit.setText(content)
        self.leStatus.setText(f"Loaded file : {fileName}")
    except Exception as e:
      self.text_edit.setText(f"Failed to load file:\n{e}")
      self.leStatus.setText(f"Failed to load file:\n{e}")
  
  def SaveFile(self):
    file_path = self.leFileName.text()
    with open(file_path, 'w') as file:
      file.write(self.text_edit.toPlainText())
      self.leStatus.setText(f"File saved to: {file_path}")

  def MakePrograms(self):
    result = subprocess.run("cd ../Cleopatra; make;cd ../PyGUIQt6", shell=True, capture_output=True, text=True)

    print("Output:", result.stdout)
    print("Error:", result.stderr)
    print("Return Code:", result.returncode)

  def CalDWBA(self):
    
    self.MakePrograms()

    if self.chkCreateInFile.isChecked :
      aMin = self.sbAngMin.value()
      aMax = self.sbAngMax.value()
      aSize = self.sbAngSize.value()





################################################## Main
if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())