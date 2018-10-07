SRC_FILES = $(wildcard ./src/*.c)
MAIN_DIR = ./main
C_FLAGS = -std=c11 -O3

founder:
	gcc $(MAIN_DIR)/founder.main.c $(SRC_FILES) -I ./inc -o ./bin/founder.exe $(C_FLAGS)

gbranch:
	gcc $(MAIN_DIR)/gain.branch.c $(SRC_FILES) -I ./inc -o ./bin/gbranch.exe $(C_FLAGS)

mapper:
	gcc $(MAIN_DIR)/mapper.main.c $(SRC_FILES) -I ./inc -o ./bin/mapper.exe $(C_FLAGS)

trajectory:
	gcc $(MAIN_DIR)/trajectory.main.c -o ./bin/trajectory.exe $(SRC_FILES) -I ./inc $(C_FLAGS)