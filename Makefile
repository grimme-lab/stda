
PROG    = stda

 FC = ifort 
# FC = lfc
 CC = gcc

### multithread ###
 LINKER = ifort -static -qopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include 
  LIBS = $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm

### sequential ###
# LINKER = ifort -static
# LIBS = ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

 CFLAGS = -O -DLINUX
 FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#################################################
OBJS=\
     stdacommon.o stringmod.o main.o pckao.o \
     header.o intpack.o velo.o  \
     onetri.o prmat.o readl.o block.o\
     stda.o stda-rw.o sutda.o sfstda.o srpapack.o intslvm.o io.o\
     linal.o readbasa.o readbasmold.o printvec.o normalize.o\
     apbtrafo.o sosor.o readxtb.o linear_response.o molden.o print_nto.o
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
	rm -f *.o *.mod $(PROG)
