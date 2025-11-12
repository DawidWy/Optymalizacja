# Kompilator i flagi
CXX = g++
CXXFLAGS = -Wpedantic -O3 -std=c++23

# Nazwa programu
TARGET = main

# Pliki źródłowe
SRCS = main.cpp \
       matrix.cpp \
       ode_solver.cpp \
       opt_alg.cpp \
       solution.cpp \
       user_funs.cpp

# Pliki obiektowe
OBJS = $(SRCS:.cpp=.o)

# Domyślny cel
all: $(TARGET)
	@$(MAKE) clean	
# Budowanie programu
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Specjalna reguła: main.cpp zależy od opt_alg.h
main.o: main.cpp opt_alg.h
	$(CXX) $(CXXFLAGS) -c main.cpp

# Reguły dla pozostałych plików (każdy ma swój header)
matrix.o: matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) -c matrix.cpp

ode_solver.o: ode_solver.cpp ode_solver.h
	$(CXX) $(CXXFLAGS) -c ode_solver.cpp

opt_alg.o: opt_alg.cpp opt_alg.h
	$(CXX) $(CXXFLAGS) -c opt_alg.cpp

solution.o: solution.cpp solution.h
	$(CXX) $(CXXFLAGS) -c solution.cpp

user_funs.o: user_funs.cpp user_funs.h
	$(CXX) $(CXXFLAGS) -c user_funs.cpp

# Czyszczenie plików tymczasowych
clean:
	rm -f $(OBJS)

# Pomoc
help:
	@echo "Dostępne cele:"
	@echo "  make         - kompiluje projekt"
	@echo "  make clean   - usuwa pliki obiektowe i program"
