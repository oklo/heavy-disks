# diskfft build: FFTW (OpenMP-threaded) + OpenMP hydro
BREW   := /opt/homebrew
OMPDIR := $(BREW)/opt/libomp
CXX    := g++
CXXFLAGS := -O3 -march=native -std=c++17 -Isrc -I$(BREW)/include -I$(OMPDIR)/include -Xpreprocessor -fopenmp
LDFLAGS  := -L$(BREW)/lib -lfftw3_omp -lfftw3 -L$(OMPDIR)/lib -lomp

all: diskdisk test_kalnajs

diskdisk: src/diskmain.cpp src/hydro.hpp src/kalnajs.hpp
	$(CXX) $(CXXFLAGS) src/diskmain.cpp -o $@ $(LDFLAGS)

test_kalnajs: src/test_kalnajs.cpp src/kalnajs.hpp
	$(CXX) $(CXXFLAGS) src/test_kalnajs.cpp -o $@ $(LDFLAGS)

clean:
	rm -f diskdisk test_kalnajs
