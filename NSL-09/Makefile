CC = g++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_091.exe : LSN_Exercise_091.o random.o chromosome.o population.o
	$(CC) random.o chromosome.o population.o LSN_Exercise_091.o -o LSN_Exercise_091.exe
LSN_Exercise_091.o :  LSN_Exercise_091.cpp LSN_Exercise_091.h
	$(CC) -c LSN_Exercise_091.cpp -o LSN_Exercise_091.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
chromosome.o : chromosome.cpp chromosome.h
	$(CC) -c chromosome.cpp -o chromosome.o $(CFLAGS)
population.o : population.cpp population.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)

clean :
	rm *.o LSN_Exercise_091.exe
