CC    =gcc
CFLAGS=-O3 -Wall
LLIB= -lm

rotation.out : main.o my_utils.o
	$(CC) $(LLIB) -o rotation.out main.o my_utils.o

main.o : main.c
	$(CC) $(CFLAGS) -c main.c

my_utils.o : my_utils.c
	$(CC) $(CFLAGS) -c my_utils.c

clean:
	@rm -rf *.o rotation.out



