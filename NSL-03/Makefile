## Define Variables
CC = g++
CFLAGS = -Wall

## Define path to subdirectories
OBJ_DIR = object_files_dir

## Define the object files needed by the excersises
OBJ03_1 = $(OBJ_DIR)/LSN_Exercise_03_1.o $(OBJ_DIR)/random.o $(OBJ_DIR)/statistical_functions.o $(OBJ_DIR)/asset.o $(OBJ_DIR)/option.o

EX03_1.exe : $(OBJ03_1)
	$(CC) $(OBJ03_1) -o EX03_1.exe
$(OBJ_DIR)/LSN_Exercise_03_1.o : LSN_Exercise_03_1.cpp
	$(CC) -c LSN_Exercise_03_1.cpp -o $(OBJ_DIR)/LSN_Exercise_03_1.o $(CFLAGS)

$(OBJ_DIR)/random.o : random.cpp random.h
	$(CC) -c random.cpp -o $(OBJ_DIR)/random.o $(CFLAGS)
$(OBJ_DIR)/statistical_functions.o : statistical_functions.cpp statistical_functions.h
	$(CC) -c statistical_functions.cpp -o $(OBJ_DIR)/statistical_functions.o $(CFLAGS)
$(OBJ_DIR)/asset.o : asset.cpp asset.h
	$(CC) -c asset.cpp -o $(OBJ_DIR)/asset.o $(CFLAGS)
$(OBJ_DIR)/option.o : option.cpp option.h
	$(CC) -c option.cpp -o $(OBJ_DIR)/option.o $(CFLAGS)

clean :
	rm $(OBJ_DIR)/*.o *.exe
