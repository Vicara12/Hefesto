COMP_OPTS = -D_GLIBCXX_DEBUG -O2 -Wall -Werror -Wno-unused-parameter -Wextra\
           -Wno-sign-compare -std=c++11

EXE_NAME = hefesto.exe

################################################################################

# get project files
ALL_CPP := $(shell ls ./src -R)


.PHONY: test