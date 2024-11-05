#!/usr/bin/python3

import os
import time
from PyQt6.QtWidgets import (
  QGridLayout, QWidget, QCheckBox
)
from PyQt6.QtCore import QUrl

from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go

class PlotWindow(QWidget):
  def __init__(self, XsecFile):
    super().__init__()

    self.setWindowTitle("DWBA Plot")
    self.setGeometry(100, 100, 800, 600)

    self.x = []
    self.data = []
    self.headers = []
    self.read_data(XsecFile)

    self.log_scale_checkbox = QCheckBox("Use Log Scale for Y-Axis")
    self.log_scale_checkbox.setChecked(True)
    self.log_scale_checkbox.stateChanged.connect(self.plot_plotly_graph)

    self.gridline_checkbox = QCheckBox("Show Gridlines")
    self.gridline_checkbox.stateChanged.connect(self.plot_plotly_graph)

    self.showMarker_checkBox = QCheckBox("Show Markers")
    self.showMarker_checkBox.stateChanged.connect(self.plot_plotly_graph)

    self.html_file = None
    self.web_view = QWebEngineView()

    layout = QGridLayout(self)
    layout.addWidget(self.showMarker_checkBox, 0, 0)
    layout.addWidget(self.log_scale_checkbox, 0, 1)
    layout.addWidget(self.gridline_checkbox, 0, 2)
    layout.addWidget(self.web_view, 1, 0, 5, 3)

    self.plot_plotly_graph()

  def read_data(self,file_path):
    self.x = []  # List for the first column
    self.data = [] # 2D list for other columns
    self.headers = []  # List to store headers

    with open(file_path, 'r') as file:
      header_found = False  # Flag to indicate if the header has been found
      for line in file:
        # Skip lines that start with '#' and empty lines
        if line.startswith('#') or not line.strip():
          continue
        
        if not header_found:
          self.headers = line.split()  # Use the split parts as headers
          header_found = True  # Set the flag to True to skip this block in future iterations
          # print(f"ELab parts found: {elab_parts}")  # Print or process this as needed
          continue
        
        # Split the line by whitespace
        parts = line.split()
        if len(parts) > 0:  # Make sure there is at least one column
          self.x.append(float(parts[0]))  # First column
          # Append the rest of the columns to data
          if len(self.data) == 0:
            # Initialize the data array with the right number of sublists
            self.data = [[] for _ in range(len(parts) - 1)]
          for i in range(len(parts) - 1):
            self.data[i].append(float(parts[i + 1]))  # Rest of the columns
            
      # print(self.headers)

  def plot_plotly_graph(self):

    if self.html_file and os.path.exists(self.html_file):
      os.remove(self.html_file)
    
    # Create a Plotly figure
    fig = go.Figure()

    if self.showMarker_checkBox.isChecked() :
      plotStyle = 'lines+markers'
    else:
      plotStyle = 'lines'

    # Add traces for each column in data against x
    for i, y in enumerate(self.data):
      fig.add_trace(go.Scatter(x=self.x, y=y, mode=plotStyle, name=self.headers[i + 1]))  # Use headers for names

    # Update layout for better presentation
    fig.update_layout(
      xaxis_title="Angle_CM [Deg]",
      yaxis_title="Xsec [mb/sr]",
      template="plotly",
      plot_bgcolor='rgba(0,0,0,0)',  # Set plot background to transparent
      paper_bgcolor='rgba(0,0,0,0)',  # Set paper background to transparent
      legend=dict(
        x=1,          # X position (1 = far right)
        y=1,          # Y position (1 = top)
        xanchor='right',  # Anchor the legend to the right
        yanchor='top',    # Anchor the legend to the top
        bgcolor='rgba(255, 255, 255, 0.5)',  # Optional: semi-transparent background for legend
        bordercolor='rgba(0, 0, 0, 0.5)',  # Optional: border color
        borderwidth=1    # Optional: border width
      ),
      yaxis=dict(
        # linecolor='black',  # Set y-axis line color to black
        type ='log' if self.log_scale_checkbox.isChecked() else 'linear',  # Toggle y-axis scale
        gridcolor='lightgray',  # Set gridline color
        gridwidth=1,  # Set gridline width (in pixels)
        showgrid = self.gridline_checkbox.isChecked()  # Toggle gridlines for y-axis
      ),
      xaxis=dict(
        # linecolor='black',  # Set x-axis line color to black
        gridcolor='lightgray',  # Set gridline color
        gridwidth=1,  # Set gridline width (in pixels)
        showgrid = self.gridline_checkbox.isChecked()  # Toggle gridlines for x-axis as well
      ),
      margin=dict(l=40, r=40, t=40, b=40),  # Set margins to reduce empty space
      # width=800,  # Optional: set fixed width for the plot
      # height=600  # Optional: set fixed height for the plot
    )

    fig.add_shape(
      # Line with reference to the plot
      type="rect",
      xref="paper",
      yref="paper",
      x0=0,
      y0=0,
      x1=1.0,
      y1=1.0,
      line=dict(
        color="black",
        width=1,
      )
    )

    # Save the plot as an HTML file in a temporary location
    timestamp = int(time.time() * 1000)  # Unique timestamp in milliseconds
    html_file = f"/tmp/plotwindow_{timestamp}.html"
    fig.write_html(html_file)
    self.html_file = html_file  # Store for cleanup
    self.web_view.setUrl(QUrl.fromLocalFile(html_file))

  def __del__(self):
    if os.path.exists(self.html_file):
       os.remove(self.html_file)
