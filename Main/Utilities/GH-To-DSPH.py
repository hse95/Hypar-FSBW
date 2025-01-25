# This code will convert the PressurePoints.txt file from the Grasshopper to the format required by the DualSPHysics
"""
Example:

from the PressurePoints.txt file:
{57.61875, 2.375, 12.375}
{57.85625, 2.125, 12.375}
{58.09375, 1.875, 12.375}
etc.

to the PressurePoints.txt file:

POINTS
57.61875  2.375 12.375
57.85625 2.125 12.375
58.09375 1.875 12.375

"""

import os

# Get the absolute path to the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Construct the absolute file path
file_path = os.path.join(current_dir, "Measure Pressure on Clusters", "PressurePoints.txt")


# open the input file in read mode
with open(file_path, "r") as file:
    # create a list to store the coordinates
    coords = []
    # read lines from the file
    lines = file.readlines()
    # for each line
    for line in lines:
        # remove the curly braces and split by comma
        line = line.replace("{","").replace("}","") .split(",")
        # remove leading/trailing whitespaces and store the coordinates
        coords.append([coord.strip() for coord in line])

# open the same file in write mode
with open(file_path, "w") as file:
    # write the header
    file.write("POINTS\n")
    # for each coordinate in the list
    for coord in coords:
        # write the coordinate to the file
        file.write(" ".join(coord) + "\n")


