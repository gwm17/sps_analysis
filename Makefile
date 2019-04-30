CC=g++
ROOTFLAGS = `root-config --cflags`
ROOTLIBS = `root-config --glibs`
CFLAGS=-g -Wall $(ROOTFLAGS)
INCLDIR=./include
SRCDIR=./src
OBJDIR=./objs
LDFLAGS= $(ROOTLIBS)
CPPFLAGS= -I$(INCLDIR)
LDLIBS= -lSpectrum
SRC=$(wildcard $(SRCDIR)/*.cpp)
OBJS=$(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXE=./analysis

.PHONY: clean

all: $(EXE)

$(EXE): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(OBJS) $(EXE)
