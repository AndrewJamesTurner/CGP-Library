# GCC complier 
CC=gcc

# GCC complier options
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -g			turns on debugging information	
	
CFLAGS= -pedantic -Wall -g 

test: src/test.c src/cgp.c src/include/cgp.h
	@$(CC) -o test src/test.c src/cgp.c src/include/cgp.h $(CFLAGS) -lm

gettingStarted: examples/gettingStarted.c src/cgp.c src/include/cgp.h
	@$(CC) -o gettingStarted examples/gettingStarted.c src/cgp.c src/include/cgp.h $(CFLAGS) -lm

createDataSet: examples/createDataSet.c src/cgp.c src/include/cgp.h
	@$(CC) -o createDataSet examples/createDataSet.c src/cgp.c src/include/cgp.h $(CFLAGS) -lm

docs: ./src/include/cgp.h ./naturaldocs/customFiles/*
	@naturaldocs -i ./src/include -i ./naturaldocs/customFiles -o FramedHTML ./docs -p ./naturaldocs

so: src/cgp.c
	@$(CC) -c -fpic src/cgp.c $(CFLAGS) 
	@gcc -shared -o libcgp.so cgp.o -lm

clean:
	@rm -f cgp.o libcgp.so cgp.dll test gettingStarted createDataSet
