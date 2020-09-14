
CC=gcc
CFLAGS= -Wall

all : psychr 

psychr : main.o svp.o psychrometrics.o
	$(CC) $(CFLAGS) -o psychr main.o svp.o psychrometrics.o -lm

main.o : main.c psychrometrics.h
	$(CC) $(CFLAGS) -c main.c

psychrometrics.o : psychrometrics.c psychrometrics.h svp.c svp.h
	$(CC) $(CFLAGS) -c psychrometrics.c -o psychrometrics.o 

svp.o : svp.c svp.h
	$(CC) $(CFLAGS) -c svp.c -o svp.o 

clean :
	rm -f psyshr *.o

