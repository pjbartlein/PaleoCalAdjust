# demo_02_adjust_1yr

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran 
FFLAGS = -O3 -ffree-line-length-160 -Wall 

# target, object and .mod files
target = demo_02_adjust_1yr
objects = pseudo_daily_interp_subs.o mp_interp_epstein_subs.o demo_02_adjust_1yr.o
modules = pseudo_daily_interp_subs.mod mp_interp_epstein_subs.mod 

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)

pseudo_daily_interp_subs.o: $(mod)pseudo_daily_interp_subs.f90
	$(FC) $(FFLAGS) -c $(mod)pseudo_daily_interp_subs.f90

mp_interp_epstein_subs.o: $(mod)mp_interp_epstein_subs.f90
	$(FC) $(FFLAGS) -c $(mod)mp_interp_epstein_subs.f90

demo_02_adjust_1yr.o: $(prog)demo_02_adjust_1yr.f90
	$(FC) $(FFLAGS) -c $(prog)demo_02_adjust_1yr.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)	
