CXX = g++
CXXFLAGS = -lm -pedantic -g
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
	rm *.dat
	rm *.png *.svg
					
run:
	./$(NAME)

plot:
	gnuplot plot_err.gnu
	gnuplot plot_integral.gnu
	gnuplot plot_map.gnu
