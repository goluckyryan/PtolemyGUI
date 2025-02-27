CC = g++
CFLAG = -g -O0

depend = ClassIsotope.h constant.h potentials.h

ALL = test

all: $(ALL)

test: inFileCreator.cpp  $(depend)
	$(CC) $(CFLAG) test.cpp -o test 


clean:
	/bin/rm -f $(ALL)