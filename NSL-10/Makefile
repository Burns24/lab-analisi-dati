# You have to change from CC to MPI or viceversa for random.o and chromosome.o
# if you want to compile the exrcise 10.1 or 10.2
CC = g++
MPI = mpic++
CFLAGS = -Wall -O3 --std=c++11

LSN_Exercise_101.exe : LSN_Exercise_101.o random.o chromosome.o
	$(CC) random.o chromosome.o LSN_Exercise_101.o -o LSN_Exercise_101.exe
LSN_Exercise_101.o :  LSN_Exercise_101.cpp LSN_Exercise_101.h
	$(CC) -c LSN_Exercise_101.cpp -o LSN_Exercise_101.o $(CFLAGS)
LSN_Exercise_102.exe : LSN_Exercise_102.o random.o chromosome.o population.o
	$(MPI) random.o chromosome.o population.o LSN_Exercise_102.o -o LSN_Exercise_102.exe
LSN_Exercise_102.o :  LSN_Exercise_102.cpp LSN_Exercise_102.h
	$(MPI) -c LSN_Exercise_102.cpp -o LSN_Exercise_102.o $(CFLAGS)

random.o : random.cpp random.h
	$(MPI) -c random.cpp -o random.o $(CFLAGS)
chromosome.o : chromosome.cpp chromosome.h
	$(MPI) -c chromosome.cpp -o chromosome.o $(CFLAGS)
population.o : population.cpp population.h
	$(MPI) -c population.cpp -o population.o $(CFLAGS)

clean :
	rm *.o *.exe
