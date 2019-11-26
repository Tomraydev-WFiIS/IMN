CXX = gcc
CXXFLAGS = -lm
OBJS = main.o
NAME = a.out

all: $(NAME)

$(NAME): $(OBJS)
	$(CXX)  -o $@  $^ $(CXXFLAGS)


main.o: main.c
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
	gnuplot plot*
