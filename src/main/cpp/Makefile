CXX = g++
CXXFLAGS = -Wall -std=c++11 -pedantic-errors -D_POSIX_C_SOURCE=199309L
OPT = -O3
DEBUG = # -g
OBJS	= PowerSum.o StirlingPowerSum.o CentralFactorialPowerSum.o EulerPowerSum.o BernoulliPowerSum.o FaulhaberPowerSum.o
SOURCE	= PowerSum.cc StirlingPowerSum.cc CentralFactorialPowerSum.cc EulerPowerSum.cc BernoulliPowerSum.cc PowerSumMain.cc FaulhaberPowerSum.cc
HEADER	= PowerSum.h StirlingPowerSum.h CentralFactorialPowerSum.h EulerPowerSum.h BernoulliPowerSum.h FaulhaberPowerSum.h
MAIN =  PowerSumMain.o
LIB = libpowersum.a
OUT	= $(LIB) PowerSum
LFLAGS	 = -lgmpxx -lgmp

all: $(OUT)

$(LIB):  $(OBJS)
	$(AR) rcs $@ $^

PowerSum: $(MAIN) $(LIB)
	$(CXX) $(CXXFLAGS) $(OPT) $(DEBUG) -o $@ $^ $(LFLAGS)

%.o: %.cc $(HEADER)
	$(CXX) -c $(CXXFLAGS) $(OPT) $(DEBUG) -o $@ $<

clean:
	rm -f $(OBJS) $(MAIN) $(OUT)
