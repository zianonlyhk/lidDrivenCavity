# **************************************************************************** #
#                                                                              #
#                                                                              #
#    Makefile                                          Personal Website        #
#                                                     ##################       #
#    By: Zian Huang <zianhuang00@gmail.com>           || room214n.com ||       #
#                                                     ##################       #
#    Created: 2023/04/10 15:22:57 by Zian Huang                                #
#                                                                              #
# **************************************************************************** #

CPP_COMPILER = g++
ADDITIONAL_LIB_PATH = -I/opt/homebrew/include/eigen3/
OPTIMISATION = -O3
COMPILE_STD = -std=c++11

OBJECTS = $(patsubst src/%.cc,obj/%.o,$(wildcard src/*.cc))

default: $(OBJECTS)
	$(CPP_COMPILER) $(COMPILE_STD) $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) $(OBJECTS) -o bin/main

obj/main.o: src/main.cc
	$(CPP_COMPILER) $(COMPILE_STD) -c $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) -o obj/main.o src/main.cc

obj/unsteady_solver.o: src/unsteady_solver.cc
	$(CPP_COMPILER) $(COMPILE_STD) -c $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) -o obj/unsteady_solver.o src/unsteady_solver.cc

clean:
	rm bin/main $(OBJECTS)