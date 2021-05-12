# demo_03_adjust_TraCE_ts

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall -ffree-line-length-160 -Wall 

# target, object and .mod files
target = demo_03_adjust_TraCE_ts
objects = spline_subs.o pseudo_daily_interp_subs.o mp_interp_harzallah_subs.o demo_03_adjust_TraCE_ts.o 
modules = spline_subs.mod pseudo_daily_interp_subs.mod mp_interp_harzallah_subs.mod 

$(target): $(objects) 
	$(FC) $(FFLAGS) -o $(target) $(objects)
	
spline_subs.o: $(mod)spline_subs.f90
	$(FC) $(FFLAGS) -c $(mod)spline_subs.f90

pseudo_daily_interp_subs.o: $(mod)pseudo_daily_interp_subs.f90
	$(FC) $(FFLAGS) -c $(mod)pseudo_daily_interp_subs.f90

mp_interp_harzallah_subs.o: $(mod)mp_interp_harzallah_subs.f90
	$(FC) $(FFLAGS) -c $(mod)mp_interp_harzallah_subs.f90

demo_03_adjust_TraCE_ts.o: $(prog)demo_03_adjust_TraCE_ts.f90
	$(FC) $(FFLAGS) -c $(prog)demo_03_adjust_TraCE_ts.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)	
