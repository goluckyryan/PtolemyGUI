#!/usr/bin/python3

import os
import datetime
import csv
import socket
import sys
import time
from functools import partial
from PyQt6.QtWidgets import QApplication, QMainWindow, QGridLayout, QPushButton, QComboBox, QWidget, QMenu, QTextEdit, QFileDialog
from PyQt6.QtCore import Qt, QPoint
from PyQt6.QtGui import QFont


class MyWindow(QMainWindow):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ptolemy GUI")
    self.setGeometry(100, 100, 1000, 700)
    self.setMinimumSize(400, 600)

    self.text_edit = QTextEdit()
    self.text_edit.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
    font = QFont("Courier New", 8)  # You can adjust the size as needed
    self.text_edit.setFont(font)

    # self.text_edit.setFixedHeight(500)
    try:
      with open("../DWBA", 'r') as file:
        content = file.read()
        self.text_edit.setText(content)
    except Exception as e:
        self.text_edit.setText(f"Failed to load file:\n{e}")

    # self.view_file_button = QPushButton("Help")
    # self.view_file_button.clicked.connect(self.open_file_viewer)

    self.cal_button = QPushButton("Calculate DWBA")
    self.cal_button.clicked.connect(self.CalDWBA)


    # Set up the layout
    layout = QGridLayout()
    layout.addWidget(self.cal_button, 0, 0)
    layout.addWidget(self.text_edit, 0, 1, 5, 5)

    # Set up the container and layout
    container = QWidget()
    container.setLayout(layout)
    self.setCentralWidget(container)


  def CalDWBA(self):
    print("Number of Row ")


if __name__ == "__main__":
  app = QApplication(sys.argv)
  window = MyWindow()
  window.show()
  sys.exit(app.exec())