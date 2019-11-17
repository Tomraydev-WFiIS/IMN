CXX=g++
OUTPUT=test
CXXFLAGS=-g -Wall
CXXFLAGS+=-pedantic
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
all: $(OUTPUT)

$(OUTPUT): $(OBJ)
	$(CXX) $? -o $@

$(OBJ): $(SRC)
	$(CXX) $(CXXFLAGS) -c $*.cpp -o $@

.PHONY: clean run copy back val class

clean:
	rm -f $(OUTPUT) $(OBJ)

run: $(OUTPUT)
	./$(OUTPUT)

copy:
	if [ -d "./copy" ]; then rm -r copy; fi
	mkdir copy
	cp $(SRC) copy
	cp *.h copy

cb:
	rm -r copy

val:
	valgrind ./$(OUTPUT)

class:
	chmod +x classgenerator.sh
	./classgenerator.sh
