CXX = g++
CXXFLAGS = -I/home/draesdom/include -g -Wall -pedantic
LDLIBS = -lsdsl -ldivsufsort -ldivsufsort64
TARGETS = build locate
SRC1 = build.cpp index/index.h index/index.cpp
SRC2 = locate.cpp index/index.h index/index.cpp

all: $(TARGETS)

build: $(SRC1)
	$(CXX) $(CXXFLAGS) $(SRC1) $(LDLIBS) -o build

locate: $(SRC2)
	$(CXX) $(CXXFLAGS) $(SRC2) $(LDLIBS) -o locate

clean:
	rm -f $(TARGETS)
