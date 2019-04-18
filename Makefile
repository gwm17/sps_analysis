CC=g++
CFLAGS=-c -Wall `root-config --cflags --glibs`
LDFLAGS=`root-config --cflags --glibs`
OBJECTS=analysis.o fit.o background.o main.o
EXECUTABLE=analysis

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -lSpectrum

%.o: %.cpp
	$(CC) -o $@ $(CFLAGS) $^ 

.PHONY: clean

clean:
	rm ./*.o ./analysis
