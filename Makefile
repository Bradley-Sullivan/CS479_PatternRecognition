# Compiler
CC = gcc

# Compiler flags
CFLAGS = -Iinclude -O2 -g

# Linker flags
LDFLAGS = -Llib -lmeschach -lm

# Output image directory
PLOT_DIR = plots/

# Data input/output directory
DATA_DIR = data

# Source directory
SRC_DIR = src

# Object files directory
OBJ_DIR = obj

# Object files to link
OBJ_LINK = lib/gnuplot_i.o

# Source files
SOURCES = $(wildcard $(SRC_DIR)/*.c)

# Object files
OBJECTS = $(SOURCES:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

# Output binary
TARGET = PA2

# Default target
all: $(TARGET)

# Link the target binary
$(TARGET): $(OBJECTS) $(OBJ_LINK)
	$(CC) $^ -o $@ $(LDFLAGS)

# Compile the object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Create the object files directory
$(OBJ_DIR):
	mkdir -p $@
	mkdir -p $(DATA_DIR) $(PLOT_DIR)

# Clean target
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

purge:
	rm -rf $(OBJ_DIR) $(TARGET)
	rm -rf $(PLOT_DIR) $(DATA_DIR)

# Prevent make from doing something with a file named clean
.PHONY: all clean
