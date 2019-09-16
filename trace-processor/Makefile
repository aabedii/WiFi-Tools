
Optimization = -O3
uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')

ifeq ($(uname_S),Darwin)	# I'm a Mac
	CC = clang++
	C_STD = c++11
  BOOST_LIBS = -lboost_system-mt
endif
ifeq ($(uname_S),Linux)		# and I'm a PC
	CC = g++
	C_STD = c++11
  BUG_WORKAROUND= -Wl,--no-as-needed 
  BOOST_LIBS = -lboost_system -lboost_iostreams 
endif

INCLUDE_DIRS = -I /usr/local/include -I /opt/local/include
LIB_DIRS = -L /usr/local/lib -L /opt/local/lib -L /usr/local/lib
THREAD_LIB = -pthread

all: trace-processor

trace-processor : trace-processor.cpp
	$(CC)  $(Optimization) $(INCLUDE_DIRS) $(LIB_DIRS) $(BOOST_LIBS) $(THREAD_LIB) -std=$(C_STD)  $(BUG_WORKAROUND)  trace-processor.cpp -o trace-processor 

clean:
	rm trace-processor

