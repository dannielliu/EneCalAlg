ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_fast6pi.o

analysis:analysis.C $(OBJS)
	$(CC) -lRooFitCore -lRooFit $^ -o $@

src/function.o : function.C	
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/bes3plotstyle.o : bes3plotstyle.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_fastpipill.o :gepep_fastpipill.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_fast4pi.o : gepep_fast4pi.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_fast6pi.o : gepep_fast6pi.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

.PHONY:clean
clean:
	rm -f analysis *.o src/*.o

cleanall:
	rm -f analysis *.o src/*.o *.eps *.pdf *.ps
