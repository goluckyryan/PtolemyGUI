#!/usr/bin/env python3
import sys

filename = sys.argv[1]

if len(sys.argv) < 2:
    print("Error: Not enough arguments provided.")
    print("Usage: ./{sys.argv[0]} filename")
    sys.exit(1)

#####################################################
import numpy as np
import re
import matplotlib.pyplot as plt

Lmax = None
sa = None
sb = None

def extract_LmaxSaSb(file_path=filename):
    global Lmax, sa, sb

    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):

            if Lmax is None and "LMAX" in line : 
                columns = line.split()
                Lmax = int(columns[6])
                print(f"LMAX = {Lmax}")
            
            if sa != None and sb is None and "2*STR" in line:
                columns = line.split()
                sb = int(float(columns[-1]))/2.
                print(f"Spin out = {sb}")

            if sa is None and "2*STR" in line:
                columns = line.split()
                sa = int(float(columns[-1]))/2.
                print(f"Spin In = {sa}")


def extract_BoundState(file_path=filename):
    r_data, y_data = [], []

    start_line = None
    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):
            if "0    R    RL,R" in line:
                if start_line is None:
                    start_line = i+1
            
            if len(r_data) > 1 and line[0] == "0" :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    if len(columns) >= 2:
                        r_data.append(float(columns[0]))  # Convert to float
                        y_data.append(float(columns[1]))

    return [r_data, y_data]

#-------------------------------------------------------
def extract_DistortedWave(file_path=filename):

    if Lmax is None or sa is None or sb is None :
        print("Please run the extract_LmaxSaSb() first")
        return

    r_list = [0.0]
    distWaveIn = np.empty((Lmax+1, int(2*(2*sa+1))), dtype=object) # going to be Lmax * (2*sa+1)
    distWaveOut = np.empty((Lmax+1, int(2*(2*sb+1))), dtype=object)  # going to be Lmax * (2*sa+1)

    for i in range(Lmax+1):
        for j in range(int(2*(2*sa+1))):
            distWaveIn[i, j] = [0.0]
        for j in range(int(2*(2*sb+1))):
            distWaveOut[i, j] = [0.0]

    start_line = None

    l_list = []

    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):

            if "0R1=" in line:
                if start_line is None:
                    start_line = i+2
                    columns = line.split()
                    r_list.append(float(columns[1]))
                    l_list = []

            if len(l_list) > Lmax:
                start_line = None

            if start_line != None:
                if start_line <= i :
                    # print(line)
                    if line[0] == "+" :
                        l = int(line[70:72])
                        for j in range(int(2*sb+1)):
                            distWaveOut[l,2*j].append(float(line[20*j+72: 20*j+82]))
                            distWaveOut[l,2*j+1].append(float(line[20*j+82:20*j+92]))
                    else:
                        l = int(line[:4])
                        l_list.append(l)
                        for j in range(int(2*sa+1)):
                            # print(float(line[20*j+4: 20*j+14]))
                            distWaveIn[l,2*j].append(float(line[20*j+4: 20*j+14]))
                            distWaveIn[l,2*j+1].append(float(line[20*j+14:20*j+24]))

    return r_list, distWaveIn, distWaveOut                    
            

#-------------------------------------------------------
def extract_ScatAmp(file_path=filename):

    if Lmax is None or sa is None or sb is None :
        print("Please run the extract_LmaxSaSb() first")
        return

    sAmpIn = []
    sAmpOut = []

    start_line = None
    
    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):
            
            if "L   REAL D1   IMAG D1   REAL D2" in line:
                if start_line is None:
                    start_line = i+1
            
            if Lmax != None and len(sAmpIn) > Lmax :
                start_line = None
                break

            if start_line != None:
                if start_line <= i :
                    columns = line.split()
                    temp = []
                    if len(columns) >= 2:
                        if line[0] == "+" : # these line are for outgoing
                            for i in range(int(2*(2*sb+1)+1)):
                                temp.append(float(columns[i+1]))
                            sAmpOut.append(temp)
                        else:                        
                            for i in range(int(2*(2*sa+1)+1)):
                                temp.append(float(columns[i]))
                            sAmpIn.append(temp)

    return sAmpIn, sAmpOut

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
def extract_Xsec(factor:float = 10, file_path=filename):
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
                        y_data.append(float(columns[1])*factor) ## factor convert fm^2 to mb

    return [x_data, y_data]

#-------------------------------------------------------
def extract_RadialMatrix(ma:str, mb:str, file_path=filename):

    if Lmax is None or sa is None or sb is None :
        print("Please run the extract_LmaxSaSb() first")
        return

    data = []

    start_line = None

    with open(file_path, "r") as file:
        for i, line in enumerate(file, start=1):

            if "0 RADIAL MATRIX ELEMENTS," in line and ma in line and mb in line:
                if start_line is None:
                    start_line = i+2
            
            if Lmax != None and len(data) > Lmax :
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
plotID = 0
plotSMID = 0

def plot_BoundState(data):
    global plotID

    x_data, y_data = data
    plt.figure(figsize=(8, 5))
    plt.plot(x_data, y_data, marker="o", linestyle="-", color="b", label="Extracted Data")
    plt.xlabel("Radius [fm]")
    plt.ylabel("Value")
    plt.title("Bound state radial function")
    plt.grid(True)

    plotID = plotID + 1
    manager = plt.get_current_fig_manager()
    manager.window.wm_geometry(f"+{850*plotID}+50")
    
    plt.show(block=False)

def plot_RadialMatrix(data, msg:str = ""):
    global plotID

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
    plt.legend()

    plotID = plotID + 1
    manager = plt.get_current_fig_manager()
    manager.window.wm_geometry(f"+{850*plotID}+50")

    plt.show(block=False)


def plot_Xsec(data, isRuth = False):
    global plotID

    x_data, y_data = data
    # plt.figure(figsize=(8, 5))
    plt.figure()
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

    plotID = plotID + 1
    manager = plt.get_current_fig_manager()
    manager.window.wm_geometry(f"+{850*plotID}+50")

    plt.show(block=False)


def plot_SMatrix(data, spin):
    global plotSMID
    
    l_data = []
    nSpin = int(2*spin+1)

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
        axes[i].set_title(f'Real and Imaginary Parts vs L for Spin {i-spin/2.:+.1f}')
        
        # Add grid lines
        axes[i].set_xlim(-1, max(l_data)+1) 
        axes[i].set_ylim(-1.1, 1.1) 
        axes[i].set_xticks(np.arange(0, max(l_data)+3, 5))  # Set x-ticks from 0 to 20 with step 5
        axes[i].grid(True)
        
        # Show the legend
        axes[i].legend()
    
    plotSMID = plotSMID + 1
    manager = plt.get_current_fig_manager()
    manager.window.wm_geometry(f"+850+{650*plotSMID}")

    # Adjust layout to prevent overlapping subplots
    plt.tight_layout()
    plt.show(block=False)


def plot_DistortWave(r_list, dw_real, dw_imag, title:str = ""):
    global plotID

    plt.figure(figsize=(8, 5))
    plt.plot(r_list, dw_real, marker="o", linestyle="-", color="b", label="Re")
    plt.plot(r_list, dw_imag, marker="x", linestyle="-", color="r", label="Im")
    plt.xlabel("Radius [fm]")
    plt.ylabel("Value")
    plt.title("Distorted Wave radial function: " +  title)
    plt.legend()
    plt.grid(True)

    plotID = plotID + 1
    manager = plt.get_current_fig_manager()
    manager.window.wm_geometry(f"+{850*plotID}+50")

    plt.show(block=False)


##############################################################################
##############################################################################

extract_LmaxSaSb() ## must be run first

bs_data = extract_BoundState()
# plot_BoundState(bs_data)

sAmpIn, sAmpOut = extract_ScatAmp()
plot_SMatrix(sAmpIn, sa)
plot_SMatrix(sAmpOut, sb)

# elXsec_data = extract_ElasticXsec()
# plot_Xsec(elXsec_data)

JA=2.5
JB=0
j=2.5
D0 = 1.55
spinFactor=(2*JB+1)/(2*JA+1)/(2*j+1)
scalingFactor=spinFactor*D0*10

print(f"spin factor : {spinFactor}")
print(f"         D0 : {D0}")
print(f"    scaling : {scalingFactor}")

xsec_data = extract_Xsec(scalingFactor)
plot_Xsec(xsec_data)
x_data, y_data = xsec_data
for i, r in enumerate(x_data):
    # if i % 5 != 0:
    #     continue
    print(f"{{{r:7.3f}, {y_data[i]:10.7f}}},")

def plot_RadialMatrix2(ma:float, mb:float, isPlot:bool=True):
    str_a = f"+{int(2*ma):2.0f}/2"
    str_b = f"+{int(2*mb):2.0f}/2"

    radmat = extract_RadialMatrix(str_a, str_b)

    if int(2*ma)%2 == 1 :
        msg_a = f"Ja = La + {int(2*ma)}/2"
    else:
        msg_a = f"Ja = La + {int(ma)}"

    if int(2*mb)%2 == 1 :
        msg_b = f"Jb = Lb + {int(2*mb)}/2"
    else:
        msg_b = f"Jb = Jb + {int(mb)}"

    if isPlot:
        plot_RadialMatrix(radmat, f": {msg_a}, {msg_b}")

    return radmat

rList, dwIn, dwOut = extract_DistortedWave()
def plot_DW(isIncoming:bool, L:int, m:float):
    if isIncoming :
        k = int(m + sa)
        plot_DistortWave(rList, dwIn[L][2*k], dwIn[L][2*k+1], f"Incoming L={L}, J=L+{m}")
    else:
        k = int(m + sb)
        plot_DistortWave(rList, dwOut[L][2*k], dwOut[L][2*k+1], f"Outgoing L={L}, J=L+{m}")

import scipy.special as sp
import scipy.interpolate as interp
def CoulombPS(L, eta):
    return np.angle(sp.gamma(1 + L + 1j * eta))

r_list, bsW = bs_data
interp_radial = interp.interp1d(r_list, bsW, kind='cubic')

def CalRadialIntgeral(L, ma, mb, isPlot:bool = True, verbose:int = 1):

    if isPlot :
        plot_DW(True, L, ma)
        plot_DW(False, L, mb)

    etaI = 0.3997
    etaO = 0.276894

    radmat = plot_RadialMatrix2(ma, mb, isPlot)
    
    prod_re, prod_im = [], []

    total = 0
    for i, r in enumerate(rList):
        ia = int(2*(ma + sa))
        ib = int(2*(mb + sb))
        dw_a = dwIn[L][ia][i] + dwIn[L][ia+1][i] * 1j
        dw_b = dwOut[L][ib][i] + dwOut[L][ib+1][i] * 1j
        bound = 0
        if rList[i] <= max(r_list) and rList[i] >= min(r_list) :
            bound = interp_radial(r)
        total += dw_a * dw_b * bound
        prod_re.append(np.real(dw_a * dw_b * bound))
        prod_im.append(np.imag(dw_a * dw_b * bound))
        if verbose >= 1:
            print(f"{i:3d} {bound:8.5f}, {rList[i]:4.1f}, {np.real(dw_a):8.5f}, {np.imag(dw_a):8.5f}, ({np.abs(dw_a):8.5f}), {np.real(dw_b):8.5f}, {np.imag(dw_b):8.5f} | {np.real(dw_a * dw_b * bound):9.6f}, {np.imag(dw_a * dw_b * bound):9.6f}")

    total = total * 0.1 * 17./16.
    phase = 1 #np.exp( 1j * (CoulombPS(L, etaI)- CoulombPS(L, etaO)) )

    if verbose >= 2:
        print("-------------------------")
        print(f" Radial integral before CoulombPS : {total}")
        print(f"                        CoulombPS : {phase}")

    total = phase * total
    print(f"                            total : [{L}, {np.real(total):8.5f}, {np.imag(total):8.6f}]")
    print(f"                           DWUCK4 : {radmat[L]}")

    plt.figure(figsize=(8, 5))
    plt.plot(rList, prod_re, marker="o", linestyle="-", color="b", label="Real")
    plt.plot(rList, prod_im, marker="x", linestyle="-", color="r", label="Imag")
    plt.xlabel("L")
    plt.ylabel("Value")
    plt.title("Product of Radial")
    plt.grid(True)
    plt.legend()
    plt.show(block=False)

    return total

# CalRadialIntgeral(3, 1, 0.5)


# radMat = plot_RadialMatrix2(-1, -0.5, True)

# for a in radMat:
#     ll, real, imag = a
#     print(f"{{{int(ll):2d}, {real:10.7f} + {imag:10.7f} I}},")

#================================================ cal Radial matrix and Plot
# radialIn = []

# for ll in range(15):
#     radialIn.append([ll, CalRadialIntgeral(ll, 1, 0.5, False)])


# l_data, real_data, imag_data = [], [], []

# for a in radialIn:
#     l_data.append(a[0])
#     real_data.append(np.real(a[1]))
#     imag_data.append(np.imag(a[1]))

# plt.figure(figsize=(8, 5))
# plt.plot(l_data, real_data, marker="o", linestyle="-", color="b", label="Real")
# plt.plot(l_data, imag_data, marker="x", linestyle="-", color="r", label="Imag")
# plt.xlabel("L")
# plt.ylabel("Value")
# plt.title("Radial Matrix (cal)")
# plt.grid(True)
# plt.show(block=False)

# plot_RadialMatrix2(1, 1/2)


input("Press Enter to exit.")