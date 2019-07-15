SRC_FILES = $(wildcard ./src/*.c)
MAIN_DIR = ./main
C_FLAGS = -std=c11 -O3

cavity-test:
	gcc $(MAIN_DIR)/cavity_test.c $(SRC_FILES) -I ./inc -o ./bin/cavtest.exe $(C_FLAGS)

kpl-mapper:
	gcc $(MAIN_DIR)/kpl.mapper.c $(SRC_FILES) -I ./inc -o ./bin/kpl-mapper.exe $(C_FLAGS)

kpr-mapper:
	gcc $(MAIN_DIR)/kpr.mapper.c $(SRC_FILES) -I ./inc -o ./bin/kpr-mapper.exe $(C_FLAGS)

new-mapper:
	gcc $(MAIN_DIR)/new.mapper.main.c $(SRC_FILES) -I ./inc -o ./bin/new-mapper.exe $(C_FLAGS)

founder:
	gcc $(MAIN_DIR)/founder.main.c $(SRC_FILES) -I ./inc -o ./bin/founder.exe $(C_FLAGS)

gbranch:
	gcc $(MAIN_DIR)/gain.branch.c $(SRC_FILES) -I ./inc -o ./bin/gbranch.exe $(C_FLAGS)

mapper:
	gcc $(MAIN_DIR)/mapper.main.c $(SRC_FILES) -I ./inc -o ./bin/mapper.exe $(C_FLAGS)

trajectory:
	gcc $(MAIN_DIR)/trajectory.main.c -o ./bin/trajectory.exe $(SRC_FILES) -I ./inc $(C_FLAGS)
