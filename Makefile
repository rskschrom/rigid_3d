# compile options
CC = g++
CCFLAGS = -O3
CLIBS = -lfmt
LPATH = -L/usr/local/lib
IPATH = -I/home/robert/research/eigen-3.4.0
OBJS = init.o quat.o physics.o geometry.o
MAIN = main.cpp

# main program compilation
main : $(MAIN) $(OBJS)
	$(CC) $(CCFLAGS) -o $@ $^ $(IPATH) $(CLIBS)

# object files
%.o : %.cpp
	$(CC) $(CCFLAGS) -c -MMD $< $(IPATH) $(CLIBS)

# clean
clean :
	rm *.o main
	
-include *.d
