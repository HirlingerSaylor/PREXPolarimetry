OPSYS = g77
PROGRAM = main
DIR_PREFIX = /w/work1701/e03004/xin/bgrate/test

ifeq ($(OPSYS),g77)
  F77 = g77
  CERN_ROOT=/apps/cernlib/i386_rhel3/2003
  FFLAGS = -c -g
  CFLAGS = -c
  LINKFLAGS = -o
  LIBS = -L$(CERN_ROOT)/lib \
         -L/usr/lib  -lm -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib
#        -L/usr/local/farm1/cernlib/99/lib  -L/usr/lib -lm
#-lf2c
endif
fort_objs =  main.o qfsrad.o fd.o fm.o fphenom.o fyukawa.o radiate.o recoil.o rom.o sig2n.o sigdel.o sigmot.o sigqfs.o sigr1.o sigr2.o sigx.o valy.o
all_objs = $(fort_objs)
#

$(PROGRAM): $(all_objs)
	rm -f $@
	$(F77) $(LINKFLAGS) $@ $(all_objs) $(LIBS)
 
#
# Get the object files from the sources
#
%.o: %.f
	rm -f $@
	$(F77) -o $@ $(FFLAGS) $<

clean: 
	rm -rf *.o *~ 


