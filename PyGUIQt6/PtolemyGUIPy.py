#!/usr/bin/python3

import os
import datetime
import csv
import socket
import sys
import time
from PyQt6.QtWidgets import QApplication, QMainWindow, QPushButton, QComboBox, QCheckBox, QLineEdit, QLabel, QVBoxLayout, QWidget, QTabWidget, QGridLayout, QMessageBox, QFileDialog, QProgressBar
from PyQt6.QtCore import Qt, QThread, QTimer, QObject, pyqtSignal
#from functools import partial
from PyQt6.QtWidgets import QApplication, QMainWindow, QTableWidget, QTableWidgetItem, QPushButton, QVBoxLayout, QWidget


class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.setGeometry(100, 100, 1000, 200)

    # Create a table with 0 row and 8 columns
    self.table = QTableWidget(0, 9)
    self.table.setHorizontalHeaderLabels(["Reaction", "gs-spin", "orbital", "spin-pi(Ex)", "Ex", "ELab [MeV/u]", "Entrance Pot.", "Exist Pot.", "Aux"])
    
    # Set up button to add new rows
    self.add_button = QPushButton("Add Row")
    self.add_button.clicked.connect(self.add_row)

    self.cal_button = QPushButton("Calculate DWBA")

    # Set up the layout
    layout = QVBoxLayout()
    layout.addWidget(self.add_button)
    layout.addWidget(self.table)
    layout.addWidget(self.cal_button)

    # Set up the container and layout
    container = QWidget()
    container.setLayout(layout)
    self.setCentralWidget(container)

    self.add_row()

  def add_row(self):
    # Add a new row to the table
    row_position = self.table.rowCount()
    self.table.insertRow(row_position)
    # Optionally populate the new row with empty items
    for column in range(5):
      self.table.setItem(row_position, column, QTableWidgetItem(""))

    combo_box1 = QComboBox()
    combo_box1.addItems(["Option 1", "Option 2", "Option 3"])
    combo_box2 = QComboBox()
    combo_box2.addItems(["Option 1", "Option 2", "Option 3"])
    self.table.setCellWidget(row_position, 6, combo_box1)
    self.table.setCellWidget(row_position, 7, combo_box2)


if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())