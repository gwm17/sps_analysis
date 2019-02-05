CC=g++
CFLAGS=-c -g -Wall `root-config --cflags`
LDFLAGS=`root-config --glibs`
SOURCES=analysis.cpp fit.cpp main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=analysis

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@
.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
clean:
	rm ./*.o ./analysis
