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
# OPTIMISATION = -O3
OPTIMISATION = -O0

OBJECTS = obj/main.o obj/unsteady_solver.o

main: $(OBJECTS)
	$(CPP_COMPILER) $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) $(OBJECTS) -o bin/main

obj/main.o: src/main.cc
	$(CPP_COMPILER) -c $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) src/main.cc -o obj/main.o

obj/unsteady_solver.o: src/unsteady_solver.cc
	$(CPP_COMPILER) -c $(ADDITIONAL_LIB_PATH) $(OPTIMISATION) src/unsteady_solver.cc -o obj/unsteady_solver.o

clean:
	rm bin/main $(OBJECTS)