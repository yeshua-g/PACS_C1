CPPFLAGS := -I../../pacs-examples/Examples/include
LDFLAGS := -L../../pacs-examples/Examples/lib
LDLIBS := -lmuparser

# Name of your executable
TARGET := gradient_descent

# Source files
SRCS := $(wildcard *.cpp)

# Object files
OBJS := $(SRCS:.cpp=.o)

# Compiler
CXX := g++ -std=c++17

# Compiler flags
CXXFLAGS := -Wall -Wextra -std=c++11

# Rule to link the object files and create the executable
$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# Rule to compile source files
%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

# Phony target to clean up generated files
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
