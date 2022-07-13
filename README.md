# Info

This is a wrapper for Ptolemy DWBA calcalation using GUI from CERN ROOT. 

The Ptolemy reference is M. H. Macfarlane and S. C. Pieper, ANL-76-11 Rev. 1, Argonne National Laboratory, 1978.

The core program (ptolemy) only works in Linux

# Usage:

in the working directory, run PtolemyGUI

Actually the program can run anywhere. It will copy the DWBA file to the present directory.

# output:

1. DWBA.in, this is the input file for the ptolmey
2. DWBA.output, this is the output file for the ptolmey
3. DWBA.Ex.txt, this is the excitation energy
4. DWBA.root, this is a root file contains the DWBA angular distribution in TGraph

# program structure

The GUI is a simple CERN ROOT GUI. 

It read the DWBA file, creates the DWBA.in with the InfileCreator.h, run the ptolemy, read the DWBA.out with the ExtractXsec.h, save the result in DWBA.Ex.txt and DWBA.root, and then use the PlotTGraphTObjArray.h to plot in root.

The optical potentials are stored in the Cleopatra/potentials.h 

