# GCC complier 
CC=gcc

# GCC complier options
	# -std=gnu99		use C99 standards
	# -pedantic		Issue all the warnings demanded by strict ISO C
	# -lm
	# -Wall			Turns on all optional warnings which are desirable for normal code.
	# -O3			turns on all optimizations
	# -I
	# -g			turns on debugging information	
	
	
CFLAGS=  -pedantic -Wall -O3 # -pg -Wextra -pg  -Wextra -ansi -pedantic -std=c99 -pedantic -Wextra -g -pg -O3 -Werror


all: src/main.c src/libCGP.c src/libCGP.h
	@$(CC) -o test src/main.c src/libCGP.c src/libCGP.h $(CFLAGS) 

so:
	@$(CC) -c -fpic src/libCGP.c $(CFLAGS) 
	@gcc -shared -o libGCP.so libCGP.o


clean:
	rm -f $(PROG)
