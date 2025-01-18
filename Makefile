fc := caf
fcflags := -O3 -I/opt/homebrew/include -fimplicit-none -fbackslash -g -fbounds-check
fclibs := -L/opt/homebrew/lib/ -L/opt/homebrew/include -lnetcdff
rm := rm -f

proj_name := model

srcdir := src
sources := $(wildcard $(srcdir)/*.f90)

objdir := build
objects = $(sources:$(srcdir)/%.f90=$(objdir)/%.o)

.PHONY: all clean

all: $(proj_name)

$(proj_name): $(objects)
	$(fc) $(fcflags) $(fclibs) $(objects) -I$(objdir) -o $(proj_name)

$(objdir)/%.o: $(srcdir)/%.f90
	$(fc) $(fcflags) -J$(objdir) -c $< -o $@

%.o: %.mod

$(objects): $(objdir)/mod_kinds.o
$(objdir)/mod_io.o: $(srcdir)/mod_io.f90 $(objdir)/mod_log.o
$(objdir)/mod_config.o: $(srcdir)/mod_config.f90 $(objdir)/mod_log.o $(objdir)/mod_util.o
$(objdir)/mod_sync.o: $(srcdir)/mod_sync.f90 $(objdir)/mod_fields.o
$(objdir)/mod_sw_dyn.o: $(srcdir)/mod_sw_dyn.f90 $(objdir)/mod_fields.o
$(objdir)/mod_model.o: $(srcdir)/mod_model.f90 $(objdir)/mod_sw_dyn.o $(objdir)/mod_sync.o $(objdir)/mod_log.o $(objdir)/mod_writer.o $(objdir)/mod_timing.o
$(objdir)/mod_fields.o: $(srcdir)/mod_fields.f90 $(objdir)/mod_tiles.o
$(objdir)/main.o: $(objects)

clean:
	$(rm) $(proj_name) build/*