# Kompilator i flagi
CXX = g++
#CXXFLAGS = -Wpedantic -O3 -std=c++23 
CXXFLAGS = -Wpedantic -O0 -fno-omit-frame-pointer -std=c++23 -g

# Nazwa programu
TARGET = main

# Pliki źródłowe
SRCS = main.cpp \
       matrix.cpp \
       ode_solver.cpp \
       opt_alg.cpp \
       solution.cpp \
       user_funs.cpp \
	   csv.cpp \
	   utils.cpp



INC_DIR = include
SRC_DIR = src
OUT_DIR = out

SRCS_FILES = $(patsubst %.cpp, $(SRC_DIR)/%.cpp, $(SRCS))
OBJS       = $(patsubst %.cpp, $(OUT_DIR)/%.o,   $(SRCS))

# Domyślny cel
all: $(TARGET)

# Budowanie programu
$(TARGET): $(OBJS)
	$(CXX) -I$(INC_DIR) $(CXXFLAGS) -o bin/$(TARGET) $(OBJS)

$(OUT_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) -I$(INC_DIR) $(CXXFLAGS) -c $< -o $@

# Utwórz katalog wyjściowy, jeśli nie istnieje
$(OUT_DIR):
	mkdir -p $(OUT_DIR)

# Czyszczenie plików tymczasowych
clean:
	rm -f $(OBJS)

# Pomoc
help:
	@echo "Dostępne cele:"
	@echo "  make         - kompiluje projekt"
	@echo "  make clean   - usuwa pliki obiektowe i program"
