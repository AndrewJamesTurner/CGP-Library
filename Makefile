# GCC complier 
CC=gcc

# GCC complier options
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -g			turns on debugging information	
	# -pg			turns on profiling information (for gprof)
	# -w			ignore warning 
	# -fpermissive
	# -std=c++11		use c++ 2011 standard
	
CFLAGS= -pedantic -Wall -O3 -fopenmp -lm

gettingStarted: examples/gettingStarted.c src/cgp.c src/cgp.h
	@$(CC) -o gettingStarted examples/gettingStarted.c src/cgp.c $(CFLAGS)  

createDataSet: examples/createDataSet.c src/cgp.c src/cgp.h
	@$(CC) -o createDataSet examples/createDataSet.c src/cgp.c $(CFLAGS)

manipulatingChromosomes: examples/manipulatingChromosomes.c src/cgp.c src/cgp.h
	@$(CC) -o manipulatingChromosomes examples/manipulatingChromosomes.c src/cgp.c $(CFLAGS)

customNodeFunction: examples/customNodeFunction.c src/cgp.c src/cgp.h
	@$(CC) -o customNodeFunction examples/customNodeFunction.c src/cgp.c $(CFLAGS) 

customFitnessFunction: examples/customFitnessFunction.c src/cgp.c src/cgp.h
	@$(CC) -o customFitnessFunction examples/customFitnessFunction.c src/cgp.c $(CFLAGS)

customSelectionScheme: examples/customSelectionScheme.c src/cgp.c src/cgp.h
	@$(CC) -o customSelectionScheme examples/customSelectionScheme.c src/cgp.c $(CFLAGS)

customReproductionScheme:  examples/customReproductionScheme.c src/cgp.c src/cgp.h
	@$(CC) -o customReproductionScheme examples/customReproductionScheme.c src/cgp.c $(CFLAGS)

averageBehaviour: examples/averageBehaviour.c src/cgp.c src/cgp.h
	@$(CC) -o averageBehaviour examples/averageBehaviour.c src/cgp.c $(CFLAGS)

neuroEvolution: examples/neuroEvolution.c src/cgp.c src/cgp.h
	@$(CC) -o neuroEvolution examples/neuroEvolution.c src/cgp.c $(CFLAGS)
	
recurrentConnections: examples/recurrentConnections.c src/cgp.c src/cgp.h
	@$(CC) -o recurrentConnections examples/recurrentConnections.c src/cgp.c $(CFLAGS)

customES: examples/customES.c src/cgp.c src/cgp.h
	@$(CC) -o customES examples/customES.c src/cgp.c $(CFLAGS)

visualization: examples/visualization.c src/cgp.c src/cgp.h
	@$(CC) -o visualization examples/visualization.c src/cgp.c $(CFLAGS)
	
printChromoEqu: examples/printChromoEqu.c src/cgp.c src/cgp.h
	@$(CC) -o printChromoEqu examples/printChromoEqu.c src/cgp.c $(CFLAGS)

multipleThreads: examples/multipleThreads.c src/cgp.c src/cgp.h
	@$(CC) -o multipleThreads examples/multipleThreads.c src/cgp.c $(CFLAGS)

so: src/cgp.c 
	@$(CC) -c -fpic src/cgp.c $(CFLAGS)
	@$(CC) -shared -o libcgp.so cgp.o -lm -fopenmp

docs: ./src/cgp.h ./naturaldocs/customFiles/*
	@naturaldocs -i ./src -i ./naturaldocs/customFiles -o HTML ./docs -p ./naturaldocs

profile: examples/averageBehaviour.c src/cgp.c src/cgp.h
	@$(CC) -o averageBehaviour examples/averageBehaviour.c src/cgp.c $(CFLAGS) -pg 
	./averageBehaviour
	gprof averageBehaviour | ./gprof2dot.py | dot -Tsvg -o profile.svg

clean:
	@rm -f cgp.o libcgp.so cgp.dll test gettingStarted createDataSet manipulatingChromosomes customNodeFunction customFitnessFunction customSelectionScheme customReproductionScheme manipluatingChromosomes averageBehaviour neuroEvolution printChromoEqu customES visualization recurrentConnections testEigen reservoirCGPANN multipleThreads *.data *.chromo *.depend *.layout *.exe *.layout *.out *.dot *.svg *.csv tmp.aux tmp.log *.pdf *.out 
	@rm -rf obj/
