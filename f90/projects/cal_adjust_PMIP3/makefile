# cal_adjust_PMIP3

# (relative) paths to source files -- needs localization
prog = ../../main_programs/
mod = ../../modules/

# netCDF and OpenMP libraries -- needs localization
NCINC = /usr/local/Cellar/netcdf/4.6.1_2/include
NCLIB = /usr/local/Cellar/netcdf/4.6.1_2/lib
LIBS = -lnetcdff -lnetcdf -lgomp

# compiler and flags
FC = gfortran
FFLAGS = -O3 -Wall

# target, object and .mod files
target = cal_adjust_PMIP3
objects = pseudo_daily_interp_subs.o calendar_effects_subs.o GISS_orbpar_subs.o \
    GISS_srevents_subs.o month_length_subs.o cmip5_netcdf_subs.o cal_adjust_PMIP3.o
modules = pseudo_daily_interp_subs.mod calendar_effects_subs.mod GISS_orbpar_subs.mod \
    GISS_srevents_subs.mod month_length_subs.mod cmip5_netcdf_subs.mod

$(target): $(objects) 
	$(FC) $(FFLAGS) -L$(NCLIB) -o $(target) $(objects) $(LIBS)

pseudo_daily_interp_subs.o: $(mod)pseudo_daily_interp_subs.f90
	$(FC) $(FFLAGS) -c $(mod)pseudo_daily_interp_subs.f90

calendar_effects_subs.o: $(mod)calendar_effects_subs.f90
	$(FC) $(FFLAGS) -c $(mod)calendar_effects_subs.f90

GISS_orbpar_subs.o: $(mod)GISS_orbpar_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_orbpar_subs.f90

GISS_srevents_subs.o: $(mod)GISS_srevents_subs.f90
	$(FC) $(FFLAGS) -c $(mod)GISS_srevents_subs.f90

month_length_subs.o: $(mod)month_length_subs.f90
	$(FC) $(FFLAGS) -c $(mod)month_length_subs.f90

cmip5_netcdf_subs.o: $(mod)cmip5_netcdf_subs.f90 
	$(FC) $(FFLAGS) -I$(NCINC) -c $(mod)cmip5_netcdf_subs.f90 

cal_adjust_PMIP3.o: $(prog)cal_adjust_PMIP3.f90
	$(FC) $(FFLAGS) -fopenmp -c $(prog)cal_adjust_PMIP3.f90

clean: 
	rm -f $(objects)
	rm -f $(modules)

clobber: clean
	rm -f $(target)