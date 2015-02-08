ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
vpath %.h ./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_fast6pi.o gepep_4k.o gepep_kk.o gepep_kpi.o gepep_kpi2.o gepep_fkkpipi.o KsAlg.o Ks0Alg.o

all: analysis AnaSinglePart
analysis:analysis.C $(OBJS)
	@echo "linking objects..."
	$(CC) -lRooFitCore -lRooFit $^ -o $@
	@echo "cleaning trash ..."
	@#-rm *.o

$(OBJS): %.o: %.C
	@echo "make object $@"
	@g++ $(ROOTINCLUDE) $(INCDIR) -c $< -o $@

#src/function.o : function.C	
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/bes3plotstyle.o : bes3plotstyle.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_fastpipill.o :gepep_fastpipill.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_fast4pi.o : gepep_fast4pi.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_fast6pi.o : gepep_fast6pi.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_4k.o : gepep_4k.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_kk.o : gepep_kk.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_kpi.o : gepep_kpi.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_kpi2.o : gepep_kpi2.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

#src/gepep_fkkpipi.o : gepep_fkkpipi.C
#	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

AnaSinglePart:AnaSinglePart.C
	@echo "compiling $@"
	$(CC) $^ -o $@


.PHONY:clean
clean:
	-rm -f analysis *.o src/*.o

cleanall:
	-rm -f analysis *.o src/*.o *.eps *.pdf *.ps
