CC := g++
CFLAGS := -std=c++17 -Wall -Wextra

BUILDDIR := build
TARGET := INDELSeek

SOURCES := main.cpp
OBJECTS := $(patsubst %.cpp,$(BUILDDIR)/%.o,$(SOURCES))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $@

$(BUILDDIR)/%.o: %.cpp
	@mkdir -p $(BUILDDIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) $(TARGET)

