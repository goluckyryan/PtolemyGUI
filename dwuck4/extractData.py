#!/usr/bin/env python3
import sys

s = float(sys.argv[1])
filename = sys.argv[2]

if len(sys.argv) < 3:
    print("Error: Not enough arguments provided.")
    print("Usage: ./{sys.argv[0]} spin filename")
    sys.exit(1)

#####################################################
import numpy as np
import re
import numpy as np
import matplotlib.pyplot as plt

def extract_BoundState(file_path=filename):
    x_data, y_data = [], []

    start_line = None
    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):
            if "0    R    RL,R" in line:
                if start_line is None:
                    start_line = i+1
            
            if len(x_data) > 1 and line[0] == "0" :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    if len(columns) >= 2:
                        x_data.append(float(columns[0]))  # Convert to float
                        y_data.append(float(columns[1]))

    return [x_data, y_data]

#-------------------------------------------------------
def extract_ScatAmp(file_path=filename):
    sAmpIn = []
    sAmpOut = []

    start_line = None
    LMAX = None
    Spin2In = None
    Spin2Out = None
    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):

            if LMAX is None and "LMAX" in line : 
                columns = line.split()
                LMAX = int(columns[6])
                print(f"LMAX = {LMAX}")
            
            if Spin2In != None and Spin2Out is None and "2*STR" in line:
                columns = line.split()
                Spin2Out = int(float(columns[-1]))
                print(f"2*Spin out = {Spin2Out}")

            if Spin2In is None and "2*STR" in line:
                columns = line.split()
                Spin2In = int(float(columns[-1]))
                print(f"2*Spin In = {Spin2In}")

            
            if "L   REAL D1   IMAG D1   REAL D2" in line:
                if start_line is None:
                    start_line = i+1
            
            if LMAX != None and len(sAmpIn) > LMAX :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    temp = []
                    if len(columns) >= 2:
                        if line[0] == "+" : # these line are for outgoing
                            for i in range(2*(Spin2Out+1)+1):
                                temp.append(float(columns[i+1]))
                            sAmpOut.append(temp)
                        else:                        
                            for i in range(2*(Spin2In+1)+1):
                                temp.append(float(columns[i]))
                            sAmpIn.append(temp)

    return Spin2In, sAmpIn, Spin2Out, sAmpOut

#-------------------------------------------------------
def extract_ElasticXsec(file_path=filename):
    x_data, y_data = [], []

    start_line = None
    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):
            if "0  Theta    Sig(1)/Coul     Sigma(1)" in line:
                if start_line is None:
                    start_line = i+1
            
            if len(x_data) > 1 and line[0] == "0" :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    if len(columns) >= 2:
                        x_data.append(float(columns[0]))  # Convert to float
                        y_data.append(float(columns[1]))

    return [x_data, y_data]


#-------------------------------------------------------
def extract_Xsec(file_path=filename):
    x_data, y_data = [], []

    start_line = None

    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):
            if "0 Theta Inelsig,fm**2" in line:
                if start_line is None:
                    start_line = i+1
            
            if len(x_data) > 1 and "0Tot-sig" in line :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    if len(columns) >= 2:
                        x_data.append(float(columns[0]))  # Convert to float
                        y_data.append(float(columns[1]))

    return [x_data, y_data]

#-------------------------------------------------------
def extract_RadialMatrix(ma:str, mb:str, file_path=filename):
    data = []

    start_line = None
    LMAX = None

    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):

            if LMAX is None and "LMAX" in line : 
                columns = line.split()
                LMAX = int(columns[6])

            if "0 RADIAL MATRIX ELEMENTS," in line and ma in line and mb in line:
                if start_line is None:
                    start_line = i+2
            
            if LMAX != None and len(data) > LMAX :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    l_data = int(line[:4])
                    real_data = float(line[4:14])
                    imag_data = float(line[14:])
                    data.append([l_data, real_data, imag_data])

    return data

#-------------------------------------------------------
def plot_BoundState(data):
    x_data, y_data = data
    plt.figure(figsize=(8, 5))
    plt.plot(x_data, y_data, marker="o", linestyle="-", color="b", label="Extracted Data")
    plt.xlabel("Radius [fm]")
    plt.ylabel("Value")
    plt.title("Bound state radial function")
    plt.grid(True)
    plt.show(block=False)

def plot_RadialMatrix(data, msg:str = ""):
    l_data, real_data, imag_data = [], [], []
    
    for a in data:
        l_data.append(a[0])
        real_data.append(a[1])
        imag_data.append(a[2])
    
    plt.figure(figsize=(8, 5))
    plt.plot(l_data, real_data, marker="o", linestyle="-", color="b", label="Real")
    plt.plot(l_data, imag_data, marker="x", linestyle="-", color="r", label="Imag")
    plt.xlabel("L")
    plt.ylabel("Value")
    plt.title("Radial Matrix " + msg)
    plt.grid(True)
    plt.show(block=False)


def plot_Xsec(data, isRuth = False):
    x_data, y_data = data
    plt.figure(figsize=(8, 5))
    plt.plot(x_data, y_data, linestyle="-", color="b", label="Extracted Data")
    plt.xlabel("Angle [deg]")
    if isRuth:
        plt.ylabel("d.s.c / Ruth ")
    else:
        plt.ylabel("d.s.c [mb/sr]")
    plt.yscale('log') 
    plt.xlim(-5, 185)
    plt.xticks(np.arange(0, 181, 20))
    plt.grid(True)
    plt.show(block=False)


def plot_SMatrix(data, spin2):
    
    l_data = []
    nSpin = spin2+1

    fig, axes = plt.subplots(1, nSpin, figsize=(6*nSpin, 4))

    smIn = []
    for a in data:
        l_data.append(int(a[0]))
        tempSM = []
        for i in range(nSpin):
            real = a[2*i+1]
            imag = a[2*i+2]
            tempSM.append((real + imag * 1j)*2j + 1 )
        smIn.append(tempSM)
        
    for i in range(nSpin):
        real_parts = [entry[i].real for entry in smIn]
        imag_parts = [entry[i].imag for entry in smIn]

        start = 1
        if i == nSpin -1 :
            start = 0

        # Plot real part vs L
        axes[i].plot(l_data[start:], real_parts[start:], label=f'Re', marker='o')
        
        # Plot imaginary part vs L
        axes[i].plot(l_data[start:], imag_parts[start:], label=f'Im', marker='x')
        
        # Adding labels and title
        axes[i].set_xlabel('L')
        axes[i].set_ylabel('Value')
        axes[i].set_title(f'Real and Imaginary Parts vs L for Spin {i-spin2/2.:+.1f}')
        
        # Add grid lines
        axes[i].set_xlim(-1, max(l_data)+1) 
        axes[i].set_ylim(-1.1, 1.1) 
        axes[i].set_xticks(np.arange(0, max(l_data)+3, 5))  # Set x-ticks from 0 to 20 with step 5
        axes[i].grid(True)
        
        # Show the legend
        axes[i].legend()
    
    # Adjust layout to prevent overlapping subplots
    plt.tight_layout()
    plt.show(block=False)

##############################################################################
##############################################################################

# bs_data = extract_BoundState()
# plot_BoundState(bs_data)

# spin2In, sAmpIn, spin2Out, sAmpOut = extract_ScatAmp()

# plot_SMatrix(sAmpIn, spin2In)
# plot_SMatrix(sAmpOut, spin2Out)

# elXsec_data = extract_ElasticXsec()
# plot_Xsec(elXsec_data)

# xsec_data = extract_Xsec()
# plot_Xsec(xsec_data)

radmat0 = extract_RadialMatrix("+ 2/2", "+ 1/2")
plot_RadialMatrix(radmat0, "+1,+1/2")

radmat1 = extract_RadialMatrix("+ 2/2", "+-1/2")
plot_RadialMatrix(radmat1, "+1,-1/2")

radmat2 = extract_RadialMatrix("+ 0/2", "+ 1/2")
plot_RadialMatrix(radmat2, "+0,+1/2")

radmat3 = extract_RadialMatrix("+ 0/2", "+-1/2")
plot_RadialMatrix(radmat3, "+0,-1/2")

radmat4 = extract_RadialMatrix("+-2/2", "+ 1/2")
plot_RadialMatrix(radmat4, "-1,+1/2")

radmat5 = extract_RadialMatrix("+-2/2", "+-1/2")
plot_RadialMatrix(radmat5, "-1,-1/2")


input("Press Enter to exit.")