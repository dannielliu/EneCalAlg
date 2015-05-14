ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include -I$(shell echo ${G4INCLUDE})
vpath %.h ./include
vpath %.o obj
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_4k.o gepep_kk.o gepep_kpi.o gepep_kpi2.o gepep_kpipi.o gepep_fkkpipi.o gepep_pipipp.o Ks0Alg.o
OBJS2= function.o gepep_fastpipill_check.o gepep_fast4pi.o gepep_fast6pi.o  gepep_kk_check.o gepep_npi_check.o gepep_kpipi_check.o gepep_fkkpipi_check.o gepep_pipipp.o Ks0Alg_check.o gepep_kpi_check.o
OBJS := $(addprefix obj/,$(OBJS))
OBJS2 := $(addprefix obj/,$(OBJS2))

all: analysis checkf AnaSinglePart AnaSomePart CompareSets SimplifyRoot MergeFiles

analysis:analysis.C $(OBJS)
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@
	@#echo "cleaning trash ..."
	@#-rm *.o

checkf:checkf.C ${OBJS2}
	@echo "compling check algorithm..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@

#$(OBJS): %.o: %.C
obj/%.o: %.C
	@echo "making object $@"
	@g++ $(ROOTINCLUDE) $(INCDIR) -c $< -o $@

%: %.C
	@echo "compiling $@"
	$(CC) $^ -o $@

#AnaSinglePart:AnaSinglePart.C
#	@echo "compiling $@"
#	$(CC) $^ -o $@

AnaSomePart:CombineParts.C
	@echo "compiling $@"
	$(CC) $^ -o $@

.PHONY:clean
clean:
	-rm -f *.o src/*.o obj/*.o analysis

cleanall:
	-rm -f analysis *.o src/*.o *.eps *.pdf *.ps
