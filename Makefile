# GCC complier 
CC=gcc

# GCC complier options
	# -std=gnu99		use C99 standards
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -g			turns on debugging information	
	
	
CFLAGS=  -pedantic -Wall -O3 


all: src/main.c src/libCGP.c src/libCGP.h
	@$(CC) -o test src/main.c src/libCGP.c src/libCGP.h $(CFLAGS) 

so:
	@$(CC) -c -fpic src/libCGP.c $(CFLAGS) 
	@gcc -shared -o libGCP.so libCGP.o


clean:
	rm -f $(PROG)
