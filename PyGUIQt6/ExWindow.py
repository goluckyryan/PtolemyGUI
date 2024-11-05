#!/usr/bin/python3

import os
import time
from PyQt6.QtWidgets import (
  QVBoxLayout, QWidget, QCheckBox
)
from PyQt6.QtCore import QUrl

from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go

from IAEANuclearData import IsotopeClass

class ExWindow(QWidget):
  def __init__(self):
    super().__init__()

    self.setWindowTitle("Ex Plot")
    self.setGeometry(100, 100, 400, 800)

    self.ASym = ""
    self.maxEx = 0
    self.data = None
    self.Iso = IsotopeClass()

    self.html_file = None
    self.web_view = QWebEngineView()

    layout = QVBoxLayout(self)
    layout.addWidget(self.web_view)

    # self.plot_Ex_graph()

  def GetEx(self, ASym :str, maxEx :float):
    self.ASym = ASym
    self.maxEx = maxEx
    self.data = self.Iso.GetExList(ASym, maxEx)    

  def plot_Ex_graph(self):

    if self.html_file and os.path.exists(self.html_file):
      os.remove(self.html_file)
    
    xShift = 0
    fontSize = 14
    plotHeight = 700
    plotWidth = 350
    yMin = -1
    
    # Create a Plotly figure
    fig = go.Figure()

    ex=self.data['energy']/1000.
    jp=self.data['jp']
    fig = go.Figure()
    fig.update_layout(plot_bgcolor='white', width=plotWidth, height = plotHeight, margin=dict(l=0, r=0, t=0, b=0))
    fig.update_layout(showlegend=False)
    fig.update_xaxes(showline=False,  visible= False, range=[-1, 3])
    fig.update_yaxes(showline=True,  visible= True, range=[yMin, self.maxEx+1])

    l=ex.last_valid_index()

    fontSizeMeV=fontSize/plotHeight*(self.maxEx+1-yMin)
    #print(fontSizeMeV)
    #adjust text label y-pos
    ypos = ex.copy()

    noOverlap = False
    loop = 0

    while noOverlap == False and loop < 2*l :
      #print("================= %d" % loop)
      for i in range(1, l+1) :
        diff = ypos[i] - ypos[i-1]
        #print("%2d | %.3f, %.3f | %.4f" % (i, ypos[i], ypos[i-1], diff))
        if diff < fontSizeMeV :
          ypos[i-1] += (diff - fontSizeMeV)/2
          ypos[i] += (fontSizeMeV - diff)/2
          if( ypos[i-1] < yMin + fontSizeMeV/2) :
            ypos[i-1] = yMin + fontSizeMeV/2
            ypos[i] = ypos[i-1] + fontSizeMeV
        #print("   | %.3f, %.3f" % (ypos[i], ypos[i-1]))

      #print(ypos)
      ###=======inspection
      count = 0
      for i in range(1, l+1) :
        diff = ypos[i] - ypos[i-1]
        if diff > fontSizeMeV :
          count = count +1

      if count == l :
        noOverlap = True

      loop += 1

    for i in range(0,l+1):
      fig.add_trace(go.Scatter(x=[xShift,1 + xShift], y=[ex[i],ex[i]],mode='lines',line=dict(color='black', width=1)))
      fig.add_trace(go.Scatter(x=[1.03 + xShift,1.1 + xShift, 1.19 + xShift], y=[ex[i],ypos[i],ypos[i]],mode='lines',line=dict(color='gray', width=1)))
      fig.add_annotation(x=1.2 + xShift, y=ypos[i], text=("%.3f, %s" % (ex[i], jp[i])), xanchor='left', font=dict(size=fontSize), showarrow=False)

    fig.add_annotation(x=0.5 + xShift, y=-0.5, text=self.ASym, font=dict(size=1.5*fontSize), showarrow=False)

    # Save the plot as an HTML file in a temporary location
    timestamp = int(time.time() * 1000)  # Unique timestamp in milliseconds
    html_file = f"/tmp/Exwindow_{timestamp}.html"
    fig.write_html(html_file)
    self.html_file = html_file  # Store for cleanup
    self.web_view.setUrl(QUrl.fromLocalFile(html_file))

  def __del__(self):
    if self.html_file and os.path.exists(self.html_file):
       os.remove(self.html_file)
