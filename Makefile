fc := caf
fcflags := -O3 -I/opt/homebrew/include
fclibs := -L/opt/homebrew/lib/ -L/opt/homebrew/include -lnetcdff
rm := rm -f

proj_name := tsunami

src := src
sources := $(wildcard $(src)/*.f90)

obj := build
objects = $(sources:$(src)/%.f90=$(obj)/%.o)

.PHONY: all clean

all: $(proj_name)

$(proj_name): $(objects)
	$(fc) $(fcflags) $(fclibs) $(objects) -I$(obj) -o $(proj_name)

$(obj)/%.o: $(src)/%.f90
	$(fc) $(fcflags) -J$(obj) -c $< -o $@

%.o: %.mod

$(objects): $(obj)/mod_kinds.o
$(obj)/mod_model.o: $(src)/mod_model.f90 $(obj)/mod_log.o
$(obj)/mod_field.o: $(src)/mod_field.f90 $(obj)/mod_io.o $(obj)/mod_parallel.o $(obj)/mod_diff.o

clean:
	$(rm) $(proj_name) build/*