
PROG    = g2molden

OSTYPE=LINUXI

ifeq ($(OSTYPE),LINUXI)
   FC = ifort 
  # FC = lfc
   CC = gcc

  LINKER = ifort -static  
 LIBS =  $(MKLROOT)/lib/intel64/libmkl_blas95_lp64.a $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
 CFLAGS = -O -DLINUX
 FFLAGS = -O3 -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include 

 
endif

ifeq ($(OSTYPE),MACOS)
    FC = ifort
    CC = gcc 
    LINKER = ifort  -qopenmp  -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
    LIBS =  ${MKLROOT}/lib/libmkl_blas95_lp64.a ${MKLROOT}/lib/libmkl_intel_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl     
    PREFLAG = -E -P
    FFLAGS = -O3 -qopenmp -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include #-check all
    FFLAGS = -O3  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include 
    CCFLAGS = -O3 -DLINUX
endif

#################################################
OBJS=\
     stringmod.o main.o 
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
