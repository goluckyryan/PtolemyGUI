CC=g++

# all: Isotope InFileCreator ExtractXSec ExtractXSecFromText PlotTGraphTObjArray 
all: InFileCreator

#Cleopatra: Cleopatra.C ../Simulation/Isotope.h ../Simulation/constant.h potentials.h InFileCreator.h ExtractXSec.h PlotTGraphTObjArray.h
#	$(CC) Cleopatra.C -o Cleopatra `root-config --cflags --glibs`

InFileCreator: InFileCreator.C InFileCreator.h ../Cleopatra/Isotope.h ../Cleopatra/constant.h potentials.h
	$(CC) InFileCreator.C -o InFileCreator `root-config --cflags --glibs`

# ExtractXSec: ExtractXSec.C ExtractXSec.h
# 	$(CC) ExtractXSec.C -o ExtractXSec `root-config --cflags --glibs`

# ExtractXSecFromText: ExtractXSecFromText.C ExtractXSec.h
# 	$(CC) ExtractXSecFromText.C -o ExtractXSecFromText `root-config --cflags --glibs`

# PlotTGraphTObjArray: PlotTGraphTObjArray.C PlotTGraphTObjArray.h
# 	$(CC) PlotTGraphTObjArray.C -o PlotTGraphTObjArray `root-config --cflags --glibs`

# Isotope: ../Cleopatra/Isotope.h ../Cleopatra/Isotope.C
# 	$(CC) Isotope.C -o Isotope
