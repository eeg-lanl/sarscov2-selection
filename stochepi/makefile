CC = g++ -std=c++17 -Isrc
OPTIM_FLAG = -g -O3
CC_FLAGS = $(OPTIM_FLAG) -Wall -c
LD_FLAGS = $(OPTIM_FLAG) -Wall
LD_LIBS = -lgsl -lgslcblas -lm -pthread
## adding -static to LD_LIBS should make the exe portable
CPP_FILES := $(wildcard src/*.cpp)
HPP_FILES := $(wildcard src/*.hpp)
CPP_SYSFILES := $(wildcard src/systems/*.cpp)
HPP_SYSFILES := $(wildcard src/systems/*.hpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
OBJ_SYSFILES := $(addprefix obj/systems/,$(notdir $(CPP_SYSFILES:.cpp=.o)))

all: obj estavoir data

clean:
	rm -rf obj
	rm -f estavoir

obj:
	mkdir -p obj
	mkdir -p obj/systems

data:
	mkdir -p data

estavoir: $(OBJ_FILES) $(OBJ_SYSFILES)
	$(CC) $(LD_FLAGS) $^ $(LD_LIBS) -o $@

obj/%.o: src/%.cpp $^
	$(CC) $(CC_FLAGS) $< -o $@

CC_FLAGS += -MMD
-include $(OBJ_FILES:.o=.d)
