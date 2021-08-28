
## Define Variables
CC = g++
CFLAGS = -Wall -O3 --std=c++11

## Define path to subdirectories
OBJ_DIR = object_files_dir

## Define the object files needed by the excersises
OBJ051 = $(OBJ_DIR)/LSN_Exercise_051.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o

EX051.exe : $(OBJ051)
	$(CC) $(OBJ051) -o EX051.exe
$(OBJ_DIR)/LSN_Exercise_051.o : LSN_Exercise_051.cpp
	$(CC) -c LSN_Exercise_051.cpp -o $(OBJ_DIR)/LSN_Exercise_051.o $(CFLAGS)

$(OBJ_DIR)/random.o : random.cpp random.h
	$(CC) -c random.cpp -o $(OBJ_DIR)/random.o $(CFLAGS)
$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)

clean :
	rm $(OBJ_DIR)/*.o *.exe
