ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_fast6pi.o gepep_kk.o gepep_kpi.o gepep_kpi2.o gepep_fkkpipi.o Pars.o

all: analysis CreateSample testsample

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

src/gepep_kk.o : gepep_kk.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_kpi.o : gepep_kpi.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_kpi2.o : gepep_kpi2.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

src/gepep_fkkpipi.o : gepep_fkkpipi.C
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

CreateSample:CreateSample.C
	$(CC) $^ -o $@

testsample:testsample.C src/function.o
	$(CC) -lRooFitCore -lRooFit $^ -o $@

src/Pars.o:Pars.cc
	g++ $(ROOTINCLUDE) $(INCDIR) -c $^ -o $@

.PHONY:clean
clean:
	rm -f analysis *.o src/*.o

cleanall:
	rm -f analysis *.o src/*.o *.eps *.pdf *.ps
