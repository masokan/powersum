CC = gcc
CCFLAGS = -Wall -std=c99 -pedantic-errors -D_POSIX_C_SOURCE=199309L
OPT = -O3
DEBUG = # -g
OBJS	= power_sum.o stirling_power_sum.o central_factorial_power_sum.o \
          euler_power_sum.o bernoulli_power_sum.o faulhaber_power_sum.o
SOURCE	= power_sum.c stirling_power_sum.c power_sum_main.c \
          euler_power_sum.c central_factorial_power_sum.c \
          bernoulli_power_sum.c faulhaber_power_sum.c
HEADER	= power_sum.h stirling_power_sum.h central_factorial_power_sum.h \
          euler_power_sum.h bernoulli_power_sum.h faulhaber_power_sum.h
AR = ar
MAIN = power_sum_main.o
LIB   = libpowersum.a
OUT	= $(LIB) powersum
LFLAGS	 = -lgmp

all: $(OUT)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

powersum: $(MAIN) $(LIB)
	$(CC) $(CCFLAGS) $(OPT) $(DEBUG) -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADER)
	$(CC) -c $(CCFLAGS) $(OPT) $(DEBUG) -o $@ $<

clean:
	rm -f $(OBJS) $(MAIN) $(OUT)
