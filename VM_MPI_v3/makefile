# choose compiler: icc or g++
CXX = mpicxx
CXXFLAGS += -O2 -std=c++11 -I${NC_DIR}/include -I${COMMON_DIR}
#CXXFLAGS += -O2 -std=c++11

LDFLAGS += -L${COMMON_DIR} -L${NC_DIR}/lib -L${H5_DIR}/lib  -Wl,-rpath=${H5_DIR}/lib
# lib path
LIBS = -lcomn -lnetcdf -lhdf5 -lhdf5_hl -lcurl -lm
# source
SRC = $(wildcard *.cpp)
# source object
SRC_OBJ = $(SRC:.cpp=.o)
# source binary
SRC_BIN =main.out

$(SRC_BIN): $(SRC_OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	rm *o *.out
