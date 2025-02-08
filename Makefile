fc := caf
fcflags := -O3 -I/opt/homebrew/include -fimplicit-none -fbackslash -g
debug_flags := -Dno_nans -fbacktrace -fbounds-check -fcheck=all -ffpe-trap=zero,invalid,overflow,underflow -Waliasing -Wcharacter-truncation -Wimplicit-interface -Wintrinsics-std -Wline-truncation -Wintrinsic-shadow -Wconversion -Wsurprising -Wunused -std=gnu -fmax-errors=5 -Warray-bounds -Wextra -Wall -pedantic
fclibs := -L/opt/homebrew/lib/ -L/opt/homebrew/include -lnetcdff
rm := rm -f

proj_name := model

srcdir := src
sources := $(wildcard $(srcdir)/*.f90)

objdir := build
objects = $(sources:$(srcdir)/%.f90=$(objdir)/%.o)

.PHONY: all clean

all: $(proj_name)
debug: fcflags += $(debug_flags)
debug: $(proj_name)

$(proj_name): $(objects)
	$(fc) $(fcflags) $(fclibs) $(objects) -I$(objdir) -o $(proj_name)

$(objdir)/%.o: $(srcdir)/%.f90
	$(fc) $(fcflags) -J$(objdir) -c $< -o $@

%.o: %.mod

$(objects): $(objdir)/mod_kinds.o
$(objdir)/mod_netcdf.o: $(srcdir)/mod_netcdf.f90 $(objdir)/mod_log.o
$(objdir)/mod_output.o: $(srcdir)/mod_output.f90 $(objdir)/mod_netcdf.o
$(objdir)/mod_config.o: $(srcdir)/mod_config.f90 $(objdir)/mod_log.o $(objdir)/mod_util.o
$(objdir)/mod_sync.o: $(srcdir)/mod_sync.f90 $(objdir)/mod_fields.o
$(objdir)/mod_sw_dyn.o: $(srcdir)/mod_sw_dyn.f90 $(objdir)/mod_fields.o $(objdir)/mod_sync.o
$(objdir)/mod_model.o: $(srcdir)/mod_model.f90 $(objdir)/mod_sw_dyn.o $(objdir)/mod_sync.o $(objdir)/mod_log.o $(objdir)/mod_input.o $(objdir)/mod_output.o $(objdir)/mod_timing.o
$(objdir)/mod_fields.o: $(srcdir)/mod_fields.f90 $(objdir)/mod_tiles.o $(objdir)/mod_netcdf.o
$(objdir)/main.o: $(objects)

clean:
	$(rm) $(proj_name) build/*