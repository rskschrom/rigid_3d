# compile options
CC = g++
CCFLAGS = 
CLIBS = 
LPATH = -L/usr/local/lib
IPATH = -I/usr/local/include
OBJS = init.o physics.o geometry.o
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
