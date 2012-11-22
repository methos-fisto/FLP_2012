OPTS  = -Wall -Wextra -g -l glpk
OPTI  = -O3 -pipe
HEADS = result.cpp flp.cpp
OBJS  = *.o
NAME  = solveFLP

launch : head
	g++ $(OBJS) $(OPTS) -o $(NAME)

head : $(HEADS)
	g++ $(HEADS) $(OPTS) -c 

opti : header
	g++  $(OBJS) $(OPTS) $(OPTI) -o $(NAME)

header : $(HEADS)
	g++ $(OPTS) $(OPTI) -c $(HEADS)

clean :
	rm -rf $(OBJS) $(NAME)

cleaner :
	rm -rf $(OBJS) $(NAME) *~ uls/*~ cls/*~ mcls/*~
