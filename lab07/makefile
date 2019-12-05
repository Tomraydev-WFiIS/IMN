CXX = g++
CXXFLAGS = -lm  -pedantic -g
OBJS = main.o
NAME = a.out

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX)  -o $@  $^ $(CXXFLAGS)


main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
			
.PHONY: clean

clean:
	rm -f *o
	rm -f ./$(NAME)
	rm *.dat || true
	rm *.png || true
	rm *.svg || true
	rm *.txt || true
					
run:
	./$(NAME)

plot:
	echo *.plot | xargs -n 1 gnuplot