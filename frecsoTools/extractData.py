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

def extract_s_matrix_data(filename):
    pattern = r"S-matrix\s+\d+\s*=\s*([\-\d\.]+)\s+([\-\d\.]+)\s+for\s+L=\s*(\d+),\s+J=\s*([\d\.]+)"

    # List to store results
    results = []

    # Variables to track the first and last "S-matrix" line
    start_line = None
    end_line = None

    # Read the file and process line by line
    with open(filename, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if "S-matrix" in line:
                # Track the first and last occurrence of "S-matrix"
                if start_line is None:
                    start_line = i
                end_line = i  # Keep updating to the latest occurrence
                
            # Extract S-matrix data if the line matches the pattern
            match = re.search(pattern, line)
            if match:
                real_part = float(match.group(1))
                imag_part = float(match.group(2))
                l_value = int(match.group(3))
                j_value = float(match.group(4))
                # Form a tuple and append to results
                results.append([real_part, imag_part, l_value, j_value])

    print(f"First 'S-matrix' at line: {start_line + 1}")
    print(f"Last 'S-matrix' at line: {end_line + 1}")

    return results

def remove_duplicates(data):
    # Sort data by L-value, then J-value
    data.sort(key=lambda x: (x[2], x[3]))  # Sort by L (index 2), then J (index 3)

    # Create a new list to store unique elements
    unique_data = []

    # Iterate over the data to remove duplicates
    for i, item in enumerate(data):
        # Add the first item or if the current item is not equal to the last added
        if i == 0 or item != data[i - 1]:
            unique_data.append(item)

    return unique_data

def group_data_by_s(data, s):
    # Initialize a dictionary for each group corresponding to J = L - s, ..., L + s
    groups = {J: [] for J in np.arange(-s, s+1, 1.0)}  # Group indices from -s to s

    
    # Iterate over the data and assign to the appropriate groups
    for entry in data:
        real_part, imaginary_part, L, J = entry
        
        # Check if J matches any value in the range L-s to L+s
        for offset in np.arange(-s, s+1, 1.0): # offset from -s to s, step 1
            if round(J, 1) == round(L + offset, 1):
                groups[offset].append(entry)

    for key in groups:
        groups[key].sort(key=lambda x: x[2])

    # Find the maximum L in each group
    max_L_in_groups = {key: max(entry[2] for entry in value) for key, value in groups.items()}

    # Find the smallest maximum L across all groups
    smallest_max_L = min(max_L_in_groups.values())

    return groups, smallest_max_L


SMdata = extract_s_matrix_data(filename)
SMdata = remove_duplicates(SMdata)

grouped_SMdata, smallest_max_L = group_data_by_s(SMdata, s)

#Mathematica input
for i in range(smallest_max_L + 1):
    print("{", end='')

    # Create a list to collect the formatted entries for each line
    line_entries = []

    for j_value, entries in grouped_SMdata.items():
        for entry in entries:
            real_part, imaginary_part, L, J = entry
            if i == int(L):
                # Format each entry and append it to the list
                line_entries.append(f"{{{L:2d}, {J:4.1f}, {real_part:9.6f}+{imaginary_part:9.6f} I}}")

    # Print the joined line without the trailing comma
    print(", ".join(line_entries), end='')

    print("},")

# Define the function to extract Ruth ratio from the file
def extract_angle_and_ratio(filename):
    angles = []
    ratios = []

    with open(filename, 'r') as file:
        # Read all lines from the file
        lines = file.readlines()

        # Find the start and end line numbers
        start_line = None
        end_line = None
        for i, line in enumerate(lines):
            if "CROSS SECTIONS FOR OUTGOING" in line:
                start_line = i + 1  # Start reading after this line
            elif "Finished all xsecs" in line:
                end_line = i
                break  # Stop once the end line is found

        # If start or end line not found, return empty lists
        if start_line is None or end_line is None:
            print("Start or end marker not found in the file.")
            return angles, ratios

        # Iterate over the identified range
        for line in lines[start_line:end_line]:
            line = line.strip()

            # Look for the lines containing "deg."
            if "deg.:" in line:
                # Extract the angle
                angle = float(line.split("deg.:")[0].strip())
                angles.append(angle)

            # Look for the lines containing "/R ="
            elif "/R =" in line:
                # Extract the /R value
                ratio = float(line.split("/R =")[1].strip())
                ratios.append(ratio)

    return angles, ratios

# Example usage
angles, ratios = extract_angle_and_ratio(filename)

# Print or save the results
print("### Differential cross section :")
for angle, ratio in zip(angles, ratios):
    if angle % 5 == 0 :
        print(f"{angle:8.1f}, {ratio:8.4f}")
print(f"{angles[-1]:8.1f}, {ratios[-1]:8.4f}")

import matplotlib.pyplot as plt


def plot_SM(groups):
    # Determine the number of groups
    num_groups = len(groups)
    
    # Create subplots (rows, columns) based on number of groups
    fig, axes = plt.subplots(1, num_groups,  figsize=(6 * num_groups,4))

    # If there's only one group, axes will not be a list, so make it iterable
    if num_groups == 1:
        axes = [axes]
    
    # Create a plot for each group
    for i, (offset, entries) in enumerate(groups.items()):
        
        # Extract real and imaginary parts, and L values
        real_parts = [entry[0] for entry in entries]
        imaginary_parts = [entry[1] for entry in entries]
        L_values = [entry[2] for entry in entries]
        
        # Plot real part vs L
        axes[i].plot(L_values, real_parts, label=f'Re', marker='o')
        
        # Plot imaginary part vs L
        axes[i].plot(L_values, imaginary_parts, label=f'Im', marker='x')
        
        # Adding labels and title
        axes[i].set_xlabel('L')
        axes[i].set_ylabel('Value')
        axes[i].set_title(f'Real and Imaginary Parts vs L for Group {offset:+.1f}')
        
        # Add grid lines
        axes[i].set_xlim(-1, max(L_values)+1) 
        axes[i].set_ylim(-1.1, 1.1) 
        axes[i].set_xticks(np.arange(0, max(L_values)+3, 5))  # Set x-ticks from 0 to 20 with step 5
        axes[i].grid(True)
        
        # Show the legend
        axes[i].legend()
    
    # Adjust layout to prevent overlapping subplots
    plt.tight_layout()

    return fig
    
def plot_Ruth(angles, ratios):
    # Plot the data
    fig = plt.figure(figsize=(8, 6))
    plt.plot(angles, ratios,  linestyle='-', color='b', label='/R values')
    plt.xscale('linear')  # Linear scale for angles
    plt.yscale('log')     # Logarithmic scale for /R values

    # Add labels, title, and legend
    plt.xlabel('CM Angle (degrees)', fontsize=12)
    plt.ylabel('d.c.s/Ruth', fontsize=12)
    plt.legend(fontsize=10)
    plt.grid(which='both', linestyle='--', linewidth=0.5)

    plt.xlim(-5, 185)
    plt.xticks(np.arange(0, 181, 20))

    # Show the plot
    plt.tight_layout()

    return fig


fig1 = plot_SM(grouped_SMdata)

fig2 = plot_Ruth(angles, ratios)

plt.show(block=False)

input("Press Enter to exit.")