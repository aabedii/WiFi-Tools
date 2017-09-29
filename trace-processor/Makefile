
Optimization = -O3
uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Darwin)	# I'm a Mac
	CC = clang++
	C_STD = c++11
endif
ifeq ($(uname_S),Linux)		# and I'm a PC
	CC = g++
	C_STD = c++11
endif

INCLUDE_DIR = /usr/local/include
LIB_DIR = /usr/local/lib
BOOST_LIBS = -lboost_system -lboost_iostreams
THREAD_LIB = -pthread
BUG_WORKAROUND= -Wl,--no-as-needed 

all: trace-processor


trace-processor : trace-processor.cpp
	$(CC)  $(Optimization) -I$(INCLUDE_DIR) -L$(LIB_DIR) $(BOOST_LIBS) $(THREAD_LIB) -std=$(C_STD)  $(BUG_WORKAROUND)  trace-processor.cpp -o trace-processor 

clean:
	rm trace-processor

