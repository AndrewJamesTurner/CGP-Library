# GCC complier 
CC=gcc

# GCC complier options
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -g			turns on debugging information	
	
CFLAGS= -pedantic -Wall 

gettingStarted: examples/gettingStarted.c src/cgp.c include/cgp.h
	@$(CC) -o gettingStarted examples/gettingStarted.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

createDataSet: examples/createDataSet.c src/cgp.c include/cgp.h
	@$(CC) -o createDataSet examples/createDataSet.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

manipulatingChromosomes: examples/manipulatingChromosomes.c src/cgp.c include/cgp.h
	@$(CC) -o manipulatingChromosomes examples/manipulatingChromosomes.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

customNodeFunction: examples/customNodeFunction.c src/cgp.c include/cgp.h
	@$(CC) -o customNodeFunction examples/customNodeFunction.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

customFitnessFunction: examples/customFitnessFunction.c src/cgp.c include/cgp.h
	@$(CC) -o customFitnessFunction examples/customFitnessFunction.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

customSelectionScheme: examples/customSelectionScheme.c src/cgp.c include/cgp.h
	@$(CC) -o customSelectionScheme examples/customSelectionScheme.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

customReproductionScheme:  examples/customReproductionScheme.c src/cgp.c include/cgp.h
	@$(CC) -o customReproductionScheme examples/customReproductionScheme.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

averageBehaviour: examples/averageBehaviour.c src/cgp.c include/cgp.h
	@$(CC) -o averageBehaviour examples/averageBehaviour.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

neuroEvolution: examples/neuroEvolution.c src/cgp.c include/cgp.h
	@$(CC) -o neuroEvolution examples/neuroEvolution.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm
	
recurrentConnections: examples/recurrentConnections.c src/cgp.c include/cgp.h
	@$(CC) -o recurrentConnections examples/recurrentConnections.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm	

customES: examples/customES.c src/cgp.c include/cgp.h
	@$(CC) -o customES examples/customES.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

visualization: examples/visualization.c src/cgp.c include/cgp.h
	@$(CC) -o visualization examples/visualization.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm
	
printChromoEqu: examples/printChromoEqu.c src/cgp.c include/cgp.h
	@$(CC) -o printChromoEqu examples/printChromoEqu.c src/cgp.c include/cgp.h $(CFLAGS) -g -lm

so: src/cgp.c 
	@$(CC) -c -fpic src/cgp.c $(CFLAGS) -O3
	@$(CC) -shared -o libcgp.so cgp.o -lm -O3

docs: ./include/cgp.h ./naturaldocs/customFiles/*
	@naturaldocs -i ./include -i ./naturaldocs/customFiles -o FramedHTML ./docs -p ./naturaldocs

clean:
	@rm -f cgp.o libcgp.so cgp.dll test gettingStarted createDataSet manipulatingChromosomes customNodeFunction customFitnessFunction customSelectionScheme customReproductionScheme manipluatingChromosomes averageBehaviour neuroEvolution printChromoEqu customES visualization *.data *.chromo *.depend *.layout *.exe *.layout *.out *.dot *.svg *.csv /obj/*
