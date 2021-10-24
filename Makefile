COMP_OPTS = -D_GLIBCXX_DEBUG -O2 -Wall -Werror -Wno-unused-parameter -Wextra\
		   -Wno-sign-compare -std=c++11

EXE_NAME = hefesto.exe

BIN_PATH     := ./bin
BUILD_PATH   := ./build
DEP_PATH     := ./dep
INCLUDE_PATH := ./include
SRC_PATH     := ./src



################################################################################

INCLUDE_PATHS := $(shell find $(INCLUDE_PATH) -type d)
ALL_CPP := $(shell find $(SRC_PATH) -type f)

.PHONY: test clean easy

all:
	g++ -o $(BIN_PATH)/$(EXE_NAME) $(ALL_CPP) -I$(INCLUDE_PATHS)

clean:
	rm $(ALL_D) $(ALL_O) $(BIN_PATH)/$(EXE_NAME)