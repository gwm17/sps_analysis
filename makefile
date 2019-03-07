CC = g++
ROOTFLAGS = `root-config --cflags --glibs`

fit_driver: fit_driver.o fit.o
	@echo "Generating fit_driver..."
	@$(CC) -o $@ $^ $(ROOTFLAGS)

%.o: %.h %.cpp
	@echo "Recompiling %.o..."
	@$(CC) -c $^ $(ROOTFLAGS)

%.o: %.cpp
	@echo "Recompiling %.o..."
	@$(CC) -c $^ $(ROOTFLAGS)

clean:
	@rm -f *~ *.o *.so *.d *.gch
