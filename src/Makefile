
CC=gcc
CFLAGS= -std=c99 -Wall

all : example chart


chart  : chart.o svp.o psychrometrics.o
	$(CC) $(CFLAGS) -o chart chart.o svp.o psychrometrics.o -lm

chart.o : chart.c psychrometrics.h
	$(CC) $(CFLAGS) -c chart.c

example : example.o svp.o psychrometrics.o
	$(CC) $(CFLAGS) -o example example.o svp.o psychrometrics.o -lm

example.o : example.c psychrometrics.h
	$(CC) $(CFLAGS) -c example.c

psychrometrics.o : psychrometrics.c psychrometrics.h svp.c svp.h
	$(CC) $(CFLAGS) -c psychrometrics.c -o psychrometrics.o 

svp.o : svp.c svp.h
	$(CC) $(CFLAGS) -c svp.c -o svp.o 

clean :
	rm -f example  *.o

