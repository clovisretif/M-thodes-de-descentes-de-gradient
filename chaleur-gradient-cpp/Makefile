CXX      = g++
CXXFLAGS = -Wall -std=c++20 -O2 -Iinclude
EXENAME  = chaleur
SRCS     = src/fonction.cpp src/main.cpp src/erreur.cpp
OBJS     = $(SRCS:src/%.cpp=build/%.o)
DEPS     = $(OBJS:.o=.d)

all: $(EXENAME)

$(EXENAME): $(OBJS)
	$(CXX) $^ -o $@

-include $(DEPS)

build/%.o: src/%.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	rm -rf build $(EXENAME)

.PHONY: all clean
