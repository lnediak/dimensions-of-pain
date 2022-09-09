INCLUDES = -I..
CXX = ccache g++ $(INCLUDES) --std=c++11 -Wall -Wextra -pedantic -g3 -O2

all: _test_convex

_test_convex: test_convex.cpp ../*.hpp
	$(CXX) -o $@ $<

clean:
	rm _test_convex

