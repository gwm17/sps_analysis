CC=g++
ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --glibs`
CFLAGS=-g -Wall $(ROOTFLAGS)
INCLDIR=./include
SRCDIR=./src
OBJDIR=./objs
PDIR=./peakfits
LDFLAGS= $(ROOTLIBS)
CPPFLAGS= -I$(INCLDIR)
LDLIBS= -lSpectrum
SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXE=./analysis
PFIT=./peakfit

.PHONY: clean all

all: $(EXE) $(PFIT)

$(EXE): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(PFIT): $(PDIR)/PeakFit.cpp
	$(CC) $(LDFLAGS) -o $@ $(LDLIBS) $(CFLAGS) $(CPPFLAGS) $^ 

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(EXE) $(PFIT)
