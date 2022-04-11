CC=clang
CFLAGS=-Wall -g

all: ruben.o mlm_init.o
	$(CC) $(CFLAGS) ruben.o mlm_init.o
	
ruben.o: ruben.c
	$(CC) $(CFLAGS) -c ruben.c
	
mlm_init.o: mlm_init.c
	$(CC) $(CFLAGS) -c mlm_init.c
	
	

clean:
	rm a.out 
	rm *.o




