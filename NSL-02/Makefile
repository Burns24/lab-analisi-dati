## Define Variables
CC = g++
CFLAGS = -Wall

## Define path to subdirectories
OBJ_DIR = object_files_dir

## Define the object files needed by the excersises
OBJ02_1 = $(OBJ_DIR)/LSN_Exercise_02_1.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o
OBJ02_2 = $(OBJ_DIR)/LSN_Exercise_02_2.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o $(OBJ_DIR)/random_walk.o

EX02_1.exe : $(OBJ02_1)
	$(CC) $(OBJ02_1) -o EX02_1.exe
$(OBJ_DIR)/LSN_Exercise_02_1.o : LSN_Exercise_02_1.cpp
	$(CC) -c LSN_Exercise_02_1.cpp -o $(OBJ_DIR)/LSN_Exercise_02_1.o $(CFLAGS)

EX02_2.exe : $(OBJ02_2)
	$(CC) $(OBJ02_2) -o EX02_2.exe
$(OBJ_DIR)/LSN_Exercise_02_2.o : LSN_Exercise_02_2.cpp
	$(CC) -c LSN_Exercise_02_2.cpp -o $(OBJ_DIR)/LSN_Exercise_02_2.o $(CFLAGS)

$(OBJ_DIR)/random.o : random.cpp random.h
	$(CC) -c random.cpp -o $(OBJ_DIR)/random.o $(CFLAGS)
$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)
$(OBJ_DIR)/random_walk.o : random_walk.cpp random_walk.h
	$(CC) -c random_walk.cpp -o $(OBJ_DIR)/random_walk.o $(CFLAGS)

clean :
	rm $(OBJ_DIR)/*.o *.exe
