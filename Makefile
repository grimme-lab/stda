PROG = std2
OSTYPE=LINUXI

#--------------------------------------------------------------------------
# see https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html details on these options

ifeq ($(OSTYPE),LINUXI)
   FC = ifx
   CC = icx

  ifdef USEILP64
  	LINKER = ifx  -Bdynamic $(CURDIR)/libcint/build/libcint.so -Bstatic
  	LIBS = -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -liomp5
   	FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/ilp64 -I$(MKLROOT)/include -i8
  else
	LINKER = ifx -Bdynamic $(CURDIR)/libcint/build/libcint.so -Bstatic
	LIBS =  -Wl,--start-group ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -liomp5
   	FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
   endif

   CFLAGS = -O -DLINUX
endif

ifeq ($(OSTYPE),MACOS)
    FC = ifx
    CC = gcc
    LINKER = ifx  -qopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
    LIBS =  ${MKLROOT}/lib/libmkl_blas95_lp64.a ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -ldl -lm
    PREFLAG = -E -P
    FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include #-check all
    FFLAGS = -O3  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include
    CCFLAGS = -O3 -DLINUX
endif

#################################################
OBJS=\
     stdacommon.o stringmod.o main.o pckao.o \
     header.o intpack.o velo.o libcint.o \
     onetri.o prmat.o readl.o block.o\
     stda.o stda-rw.o stda-rw_dual.o sutda.o sfstda.o srpapack.o intslvm.o io.o\
     linal.o readbasa.o readbasmold.o printvec.o normalize.o 2PA.o\
     apbtrafo.o sosor.o readxtb.o linear_response.o molden.o print_nto.o xstd.o full.o
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
