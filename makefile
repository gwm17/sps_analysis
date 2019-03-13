CC = g++
ROOTFLAGS = `root-config --cflags --glibs`

rootfile: fit_driver fit_driver.cpp fit.cpp fit.h
	@./$<

fit_driver: fit_driver.o fit.o
	@echo "Generating fit_driver..."
	@$(CC) -g -o $@ $^ $(ROOTFLAGS)

%.o: %.h %.cpp
	@echo "Recompiling %.o..."
	@$(CC) -g -c $^ $(ROOTFLAGS)

%.o: %.cpp
	@echo "Recompiling %.o..."
	@$(CC) -g -c $^ $(ROOTFLAGS)

clean:
	@rm -f *~ *.o *.so *.d *.gch
