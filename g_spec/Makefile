
PROG    = g_spec

 FC = ifort 
 CC = gcc
 LINKER = ifort -static  
 LIBS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
 CFLAGS = -O -DLINUX
 FFLAGS = -O3 -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include 

#################################################
OBJS=\
     g_spec.o 
#################################################

%.o: %.f90
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f
	@echo "making $@ from $<"
	$(FC) $(FFLAGS) -c $< -o $@

$(PROG):     $(OBJS) 
		@echo  "Loading $(PROG) ... "
		@$(LINKER) $(OBJS) $(LIBS) -o $(PROG)

clean:
	rm -f *.o $(PROG)
