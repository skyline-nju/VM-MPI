# choose compiler: icc or g++
CXX = mpicxx
CXXFLAGS += -O3 -std=c++11
#CXXFLAGS += -O2 -std=c++11

# LDFLAGS += -L${COMMON_DIR}
# lib path
LIBS = -lm
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN = a.out

$(SRC_BIN): $(SRC_OBJ)
	$(CXX) -o $@ $^ $(LIBS)
#	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm *o *.out
