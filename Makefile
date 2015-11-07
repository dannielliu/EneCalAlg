ROOTINCLUDE = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)
INCDIR = -I./include -I$(shell echo ${G4INCLUDE})
vpath %.h ./include
vpath %.o obj
VPATH = src
CC = g++ $(ROOTINCLUDE) $(ROOTLIB) $(INCDIR)
OBJS= function.o bes3plotstyle.o gepep_fastpipill.o gepep_fast4pi.o gepep_4k.o gepep_kk.o gepep_kpi.o gepep_ppi.o gepep_kpi2.o gepep_kpipi.o gepep_fkkpipi.o gepep_pipipp.o Ks0Alg.o mumu.o
OBJS2= function.o gepep_fastpipill_check.o gepep_fast4pi_check.o gepep_fast6pi.o  gepep_kk_check.o gepep_npi_check.o gepep_kpipi_check.o gepep_fkkpipi_check.o gepep_pipipp.o Ks0Alg_check.o gepep_kpi_check.o mumu_check.o gepep_lambdac.o
OBJS := $(addprefix obj/,$(OBJS))
OBJS2 := $(addprefix obj/,$(OBJS2))

targets0 = analysis checkf AnaSinglePart AnaSomePart CompareSets SimplifyRoot MergeFiles
targets := $(addprefix bin/,$(targets0))

all: $(targets)

bin/analysis:analysis.C $(OBJS)
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@

bin/checkf:checkf.C ${OBJS2}
	@echo "compling check algorithm..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@


k3pi: bin/analysis_k3pi

bin/analysis_k3pi: analysis_k3pi.C obj/gepep_k3pi.o obj/function.o
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@

kstar: bin/analysis_kstar
bin/analysis_kstar: analysis_kstar.C obj/gepep_kpi_Kstar.o obj/function.o
	@echo "compling analysis algorithm, linking objects..."
	$(CC) -lRooFitCore -lRooFit -lMathMore $^ -o $@

#$(OBJS): %.o: %.C
obj/%.o: %.C
	@echo "making object $@"
	@g++ $(ROOTINCLUDE) $(INCDIR) -c $< -o $@

bin/%: %.C
	@echo "compiling $@"
	$(CC) $^ -o $@

bin/AnaSomePart:CombineParts.C
	@echo "compiling $@"
	$(CC) $^ -o $@

.PHONY:clean
clean:
	-rm -f  $(OBJS) $(OBJS2) $(targets)

cleanall:
	-rm -f analysis *.o src/*.o *.eps *.pdf *.ps




