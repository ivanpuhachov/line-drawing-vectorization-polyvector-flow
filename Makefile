#let's try if that works
CC=g++-4.8
CPPFLAGS=-std=c++11 -O3 -g -fopenmp
OBJDIR=objects
SOURCES=$(wildcard *.cpp)
OBJECTS=$(patsubst %.cpp, $(OBJDIR)/%.o, $(SOURCES))
INCLUDES=-I/home/bmpix/eigen-3.3.4/ -I/opt/IpOpt/build/include/coin/ -I/opt/gurobi811/linux64/include/
LDFLAGS=`pkg-config --libs opencv` -lpthread -lcuda -Wl,-rpath-link,/usr/local/cuda-8.0/lib64/ -L/opt/IpOpt/build/lib/ -L/home/bmpix/Ipopt-3.12.4/build/lib/ -L/opt/gurobi811/linux64/lib/ -lipopt -lcoinlapack -lcoinblas -lgurobi_c++ -lgurobi_g++5.2 -lgurobi81

all: $(OBJDIR)/$(OBJECTS)
	$(CC) $(OBJECTS) $(CPPFLAGS) $(LDFLAGS) -o vect

$(OBJDIR)/%.o: ./%.cpp
	$(CC) $(CPPFLAGS) $(INCLUDES) -c $< -o $@
print-%  : ; @echo $* = $($*)
clean:
	rm -f $(OBJDIR)/*.o
	rm -f vect
