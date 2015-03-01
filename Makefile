ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include
vpath %.h ./include
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_fast6pi.o gepep_4k.o gepep_kk.o gepep_kpi.o gepep_kpi2.o gepep_fkkpipi.o Ks0Alg.o
OBJS2= function.o gepep_fastpipill_check.o gepep_fast4pi.o gepep_fast6pi.o gepep_fkkpipi_check.o Ks0Alg_check.o

all: analysis checkf AnaSinglePart AnaSomePart CompareSets SimplifyRoot MergeFiles

analysis:analysis.C $(OBJS)
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit $^ -o $@
	@#echo "cleaning trash ..."
	@#-rm *.o

checkf:checkf.C ${OBJS2}
	@echo "compling check algorithm..."
	$(CC) -lRooFitCore -lRooFit $^ -o $@

#$(OBJS): %.o: %.C
%.o: %.C
	@echo "making object $@"
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

AnaSomePart:CombineParts.C
	@echo "compiling $@"
	$(CC) $^ -o $@

CompareSets:CompareSets.C
	@echo "compiling $@"
	$(CC) $^ -o $@

SimplifyRoot:SimplifyRoot.C
	@echo "compiling $@"
	$(CC) $^ -o $@

MergeFiles:MergeFiles.C
	@echo "compiling $@"
	$(CC) $^ -o $@


.PHONY:clean
clean:
	-rm -f *.o src/*.o

cleanall:
	-rm -f analysis *.o src/*.o *.eps *.pdf *.ps
