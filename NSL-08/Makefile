CC = g++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_082.exe : LSN_Exercise_082.o random.o
	$(CC) random.o LSN_Exercise_082.o -o LSN_Exercise_082.exe
LSN_Exercise_082.o : LSN_Exercise_082.cpp
	$(CC) -c LSN_Exercise_082.cpp -o LSN_Exercise_082.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o LSN_Exercise_082.exe
