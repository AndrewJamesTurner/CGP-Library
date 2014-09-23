# GCC complier 
CC=gcc

# GCC complier options
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -g			turns on debugging information	
	# -pg			turns on profiling information (for gprof)
	
CFLAGS= -pedantic -Wall -O3

gettingStarted: examples/gettingStarted.c src/cgp.c src/cgp.h
	@$(CC) -o gettingStarted examples/gettingStarted.c src/cgp.c src/cgp.h $(CFLAGS) -lm

createDataSet: examples/createDataSet.c src/cgp.c src/cgp.h
	@$(CC) -o createDataSet examples/createDataSet.c src/cgp.c src/cgp.h $(CFLAGS)  -lm

manipulatingChromosomes: examples/manipulatingChromosomes.c src/cgp.c src/cgp.h
	@$(CC) -o manipulatingChromosomes examples/manipulatingChromosomes.c src/cgp.c src/cgp.h $(CFLAGS) -lm

customNodeFunction: examples/customNodeFunction.c src/cgp.c src/cgp.h
	@$(CC) -o customNodeFunction examples/customNodeFunction.c src/cgp.c src/cgp.h $(CFLAGS) -lm

customFitnessFunction: examples/customFitnessFunction.c src/cgp.c src/cgp.h
	@$(CC) -o customFitnessFunction examples/customFitnessFunction.c src/cgp.c src/cgp.h $(CFLAGS) -lm

customSelectionScheme: examples/customSelectionScheme.c src/cgp.c src/cgp.h
	@$(CC) -o customSelectionScheme examples/customSelectionScheme.c src/cgp.c src/cgp.h $(CFLAGS) -lm

customReproductionScheme:  examples/customReproductionScheme.c src/cgp.c src/cgp.h
	@$(CC) -o customReproductionScheme examples/customReproductionScheme.c src/cgp.c src/cgp.h $(CFLAGS) -lm

averageBehaviour: examples/averageBehaviour.c src/cgp.c src/cgp.h
	@$(CC) -o averageBehaviour examples/averageBehaviour.c src/cgp.c src/cgp.h $(CFLAGS) -lm

neuroEvolution: examples/neuroEvolution.c src/cgp.c src/cgp.h
	@$(CC) -o neuroEvolution examples/neuroEvolution.c src/cgp.c src/cgp.h $(CFLAGS) -lm
	
recurrentConnections: examples/recurrentConnections.c src/cgp.c src/cgp.h
	@$(CC) -o recurrentConnections examples/recurrentConnections.c src/cgp.c src/cgp.h $(CFLAGS) -lm	

customES: examples/customES.c src/cgp.c src/cgp.h
	@$(CC) -o customES examples/customES.c src/cgp.c src/cgp.h $(CFLAGS) -lm

visualization: examples/visualization.c src/cgp.c src/cgp.h
	@$(CC) -o visualization examples/visualization.c src/cgp.c src/cgp.h $(CFLAGS) -lm
	
printChromoEqu: examples/printChromoEqu.c src/cgp.c src/cgp.h
	@$(CC) -o printChromoEqu examples/printChromoEqu.c src/cgp.c src/cgp.h $(CFLAGS) -lm

so: src/cgp.c 
	@$(CC) -c -fpic src/cgp.c $(CFLAGS) -O3
	@$(CC) -shared -o libcgp.so cgp.o -lm -O3

docs: ./src/cgp.h ./naturaldocs/customFiles/*
	@naturaldocs -i ./src -i ./naturaldocs/customFiles -o FramedHTML ./docs -p ./naturaldocs

profile: examples/averageBehaviour.c src/cgp.c src/cgp.h
	@$(CC) -o averageBehaviour examples/averageBehaviour.c src/cgp.c src/cgp.h $(CFLAGS) -pg -lm
	./averageBehaviour
	gprof averageBehaviour | ./gprof2dot.py | dot -Tsvg -o profile.svg

python: src/cgp.c src/cgp.h
	@swig -python bindings/cgp.i
	@$(CC) -c src/cgp.c bindings/cgp_wrap.c -I/usr/include/python2.7
	@ld -shared cgp.o cgp_wrap.o -o bindings/_cgp.so
	@rm cgp.o cgp_wrap.o 
	#@rm bindings/cgp_wrap.c
clean:
	@rm -f cgp.o libcgp.so cgp.dll test gettingStarted createDataSet manipulatingChromosomes customNodeFunction customFitnessFunction customSelectionScheme customReproductionScheme manipluatingChromosomes averageBehaviour neuroEvolution printChromoEqu customES visualization recurrentConnections *.data *.chromo *.depend *.layout *.exe *.layout *.out *.dot *.svg *.csv /obj/* tmp.aux tmp.log tmp.pdf *.out
