# Makefile for s2 library
# Jason McEwen


# ======== OPTIONS ========

USEPGPLOT = no
#USEPGPLOT = yes


# ======== COMPILER ========

FC      = gfortran
#FC      = f95
#FC      = g95

ifneq ($(USEPGPLOT),yes)
  OPTPGPLOT = -DNO_PGPLOT
endif

ifeq ($(FC),f95)
  OPTF95 = -w=x95
endif

OPT = $(OPTPGPLOT) $(OPTF95) -DWMAP5 -m64
#OPT = $(OPTPGPLOT) $(OPTF95) -DMAKE_COADDED_DATA_MAP -DWMAP3 -DDEBUG


# ======== LINKS ========

PROGDIR = /Users/mcewen/src

HPIXDIR = $(PROGDIR)/Healpix
HPIXLIB = $(HPIXDIR)/lib
HPIXLIBNM= healpix
HPIXINC = $(HPIXDIR)/include

S2DIR  = $(PROGDIR)/s2
S2LIB  = $(S2DIR)/lib
S2LIBNM= s2
S2INC  = $(S2DIR)/include
S2SRC  = $(S2DIR)/src/mod
S2PROG = $(S2DIR)/src/prog
S2BIN  = $(S2DIR)/bin
S2DOC  = $(S2DIR)/doc

CFITSIOLIB   = $(PROGDIR)/cfitsio/lib
CFITSIOLIBNM = cfitsio

PGPLOTLIB    = $(PROGDIR)/pgplot
PGPLOTLIBNM  = pgplot
X11LIB       = /usr/X11R6/lib
X11LIBNM     = X11


# ======== FFFLAGS ========

FFLAGS  = -I$(HPIXINC) -I$(S2INC)


# ======== LDFLAGS ========

ifeq ($(USEPGPLOT),yes)
  LDFLAGSPGPLOT = -L$(PGPLOTLIB) -L$(X11LIB) \
                  -l$(PGPLOTLIBNM) -l$(X11LIBNM)
endif

LDFLAGS =  -L$(S2LIB) -l$(S2LIBNM) \
           -L$(HPIXLIB) -l$(HPIXLIBNM) \
           -L$(CFITSIOLIB) -l$(CFITSIOLIBNM) \
           $(LDFLAGSPGPLOT)


# ======== PPFLAGS ========

ifeq ($(FC),f95)
  PPFLAGS = -fpp $(OPT)
else ifeq ($(FC),g95)
  PPFLAGS = -cpp $(OPT)
else ifeq ($(FC),gfortran)
  PPFLAGS = -x f95-cpp-input $(OPT)
endif


# ======== OBJECT FILES TO MAKE ========

S2OBJ  = $(S2INC)/s2_types_mod.o   \
          $(S2INC)/s2_error_mod.o  \
          $(S2INC)/s2_vect_mod.o   \
          $(S2INC)/s2_distn_mod.o  \
          $(S2INC)/s2_wnoise_mod.o \
          $(S2INC)/s2_pl_mod.o     \
          $(S2INC)/s2_cmb_mod.o    \
          $(S2INC)/s2_dl_mod.o     \
          $(S2INC)/s2_sky_mod.o    \
          $(S2INC)/s2_ylm_mod.o 


# ======== MAKE RULES ========

default: all

all:     lib prog

lib:	 $(S2LIB)/lib$(S2LIBNM).a

prog:    $(S2BIN)/s2_nonzero     \
         $(S2BIN)/s2_sky2cl      \
         $(S2BIN)/s2_sky2map     \
         $(S2BIN)/s2_sky2alm     \
         $(S2BIN)/s2_map2sky     \
         $(S2BIN)/s2_map2alm     \
         $(S2BIN)/s2_map2matmap  \
         $(S2BIN)/s2_alm2sky     \
         $(S2BIN)/s2_alm2map     \
         $(S2BIN)/s2_alm2matalm  \
         $(S2BIN)/s2_plplot      \
         $(S2BIN)/s2_sky2plplot  \
         $(S2BIN)/s2_cl2ascii    \
         $(S2BIN)/s2_skyadd      \
         $(S2BIN)/s2_skydown     \
         $(S2BIN)/s2_skymultiply \
         $(S2BIN)/s2_skyoffset \
         $(S2BIN)/s2_gcmb        \
         $(S2BIN)/s2_skyup       \
         $(S2BIN)/s2_skyprod     \
         $(S2BIN)/s2_skymask     \
         $(S2BIN)/s2_grf2        \
         $(S2BIN)/s2_gcmbcoad    \
         $(S2BIN)/s2_skyrot      \
         $(S2BIN)/s2_almrot      \
         $(S2BIN)/s2_skyorder    \
         $(S2BIN)/s2_dlwrite     \
         $(S2BIN)/s2_skyerror    \
         $(S2BIN)/s2_skyrms      \
         $(S2BIN)/s2_skythres    \
         $(S2BIN)/s2_ylm         \
         $(S2BIN)/s2_nmask       \
         $(S2BIN)/s2_sky2fsht

$(S2INC)/%.o: $(S2SRC)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 
	mv *.mod $(S2INC)

$(S2INC)/%.o: $(S2PROG)/%.f90
	$(FC) $(FFLAGS) $(PPFLAGS) -c $< -o $@ 

# Library

$(S2LIB)/lib$(S2LIBNM).a: $(S2OBJ)
	ar -r $(S2LIB)/lib$(S2LIBNM).a $(S2OBJ)


# Documentation

docs:
	./f90doc_fpp $(S2SRC)/*.f90
	./f90doc_fpp $(S2PROG)/*.f90
	mv *.html $(S2DOC)/.
	./addstyle $(S2DOC)/s2_*

cleandocs:
	rm -f $(S2DOC)/s2_*.html


# Cleaning up

clean:	tidy
	rm -f $(S2INC)/*.mod
	rm -f $(S2INC)/*.o
	rm -f $(S2LIB)/lib$(S2LIBNM).a
	rm -f $(S2BIN)/*

tidy:
	rm -f *.mod
	rm -f $(S2SRC)/*~ 
	rm -f $(S2PROG)/*~ 

# Module dependencies

$(S2INC)/s2_types_mod.o: $(S2SRC)/s2_types_mod.f90
$(S2INC)/s2_error_mod.o: $(S2SRC)/s2_error_mod.f90  \
                             $(S2INC)/s2_types_mod.o
$(S2INC)/s2_sky_mod.o:   $(S2SRC)/s2_sky_mod.f90    \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o  \
                             $(S2INC)/s2_vect_mod.o   \
                             $(S2INC)/s2_dl_mod.o   \
                             $(S2INC)/s2_pl_mod.o
$(S2INC)/s2_distn_mod.o: $(S2SRC)/s2_distn_mod.f90  \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o
$(S2INC)/s2_pl_mod.o:    $(S2SRC)/s2_pl_mod.f90     \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o
$(S2INC)/s2_wnoise_mod.o:$(S2SRC)/s2_wnoise_mod.f90 \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o  \
                             $(S2INC)/s2_distn_mod.o  \
                             $(S2INC)/s2_sky_mod.o
$(S2INC)/s2_cmb_mod.o:   $(S2SRC)/s2_cmb_mod.f90    \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o  \
                             $(S2INC)/s2_distn_mod.o  \
                             $(S2INC)/s2_pl_mod.o     \
                             $(S2INC)/s2_sky_mod.o    \
                             $(S2INC)/s2_wnoise_mod.o
$(S2INC)/s2_vect_mod.o:  $(S2SRC)/s2_vect_mod.f90   \
                             $(S2INC)/s2_types_mod.o  \
                             $(S2INC)/s2_error_mod.o
$(S2INC)/s2_dl_mod.o:  $(S2SRC)/s2_dl_mod.f90   \
                             $(S2INC)/s2_types_mod.o
$(S2INC)/s2_ylm_mod.o:  $(S2SRC)/s2_ylm_mod.f90   \
                             $(S2INC)/s2_types_mod.o \
                             $(S2INC)/s2_error_mod.o \
                             $(S2INC)/s2_dl_mod.o \
                             $(S2INC)/s2_sky_mod.o    

# Program dependencies and compilation

$(S2INC)/s2_nonzero.o:     $(S2PROG)/s2_nonzero.f90 lib
$(S2BIN)/s2_nonzero:       $(S2INC)/s2_nonzero.o
	$(FC) -o $(S2BIN)/s2_nonzero $(S2INC)/s2_nonzero.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_sky2cl.o:	     $(S2PROG)/s2_sky2cl.f90 lib
$(S2BIN)/s2_sky2cl:        $(S2INC)/s2_sky2cl.o
	$(FC) -o $(S2BIN)/s2_sky2cl $(S2INC)/s2_sky2cl.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_sky2map.o:     $(S2PROG)/s2_sky2map.f90 lib
$(S2BIN)/s2_sky2map:	     $(S2INC)/s2_sky2map.o
	$(FC) -o $(S2BIN)/s2_sky2map $(S2INC)/s2_sky2map.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_sky2alm.o:     $(S2PROG)/s2_sky2alm.f90 lib
$(S2BIN)/s2_sky2alm:	     $(S2INC)/s2_sky2alm.o
	$(FC) -o $(S2BIN)/s2_sky2alm $(S2INC)/s2_sky2alm.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_map2sky.o:     $(S2PROG)/s2_map2sky.f90 lib
$(S2BIN)/s2_map2sky:	     $(S2INC)/s2_map2sky.o
	$(FC) -o $(S2BIN)/s2_map2sky $(S2INC)/s2_map2sky.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_map2alm.o:     $(S2PROG)/s2_map2alm.f90 lib
$(S2BIN)/s2_map2alm:	     $(S2INC)/s2_map2alm.o
	$(FC) -o $(S2BIN)/s2_map2alm $(S2INC)/s2_map2alm.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_alm2sky.o:     $(S2PROG)/s2_alm2sky.f90 lib
$(S2BIN)/s2_alm2sky:	     $(S2INC)/s2_alm2sky.o
	$(FC) -o $(S2BIN)/s2_alm2sky $(S2INC)/s2_alm2sky.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_alm2map.o:     $(S2PROG)/s2_alm2map.f90 lib
$(S2BIN)/s2_alm2map:	     $(S2INC)/s2_alm2map.o
	$(FC) -o $(S2BIN)/s2_alm2map $(S2INC)/s2_alm2map.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_plplot.o:	     $(S2PROG)/s2_plplot.f90 lib
$(S2BIN)/s2_plplot:	     $(S2INC)/s2_plplot.o
	$(FC) -o $(S2BIN)/s2_plplot $(S2INC)/s2_plplot.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_sky2plplot.o:  $(S2PROG)/s2_sky2plplot.f90 lib
$(S2BIN)/s2_sky2plplot:    $(S2INC)/s2_sky2plplot.o
	$(FC) -o $(S2BIN)/s2_sky2plplot $(S2INC)/s2_sky2plplot.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_cl2ascii.o:    $(S2PROG)/s2_cl2ascii.f90 lib
$(S2BIN)/s2_cl2ascii:      $(S2INC)/s2_cl2ascii.o
	$(FC) -o $(S2BIN)/s2_cl2ascii $(S2INC)/s2_cl2ascii.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyadd.o:      $(S2PROG)/s2_skyadd.f90 lib
$(S2BIN)/s2_skyadd:        $(S2INC)/s2_skyadd.o
	$(FC) -o $(S2BIN)/s2_skyadd $(S2INC)/s2_skyadd.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyprod.o:     $(S2PROG)/s2_skyprod.f90 lib
$(S2BIN)/s2_skyprod:       $(S2INC)/s2_skyprod.o
	$(FC) -o $(S2BIN)/s2_skyprod $(S2INC)/s2_skyprod.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skymask.o:     $(S2PROG)/s2_skymask.f90 lib
$(S2BIN)/s2_skymask:       $(S2INC)/s2_skymask.o
	$(FC) -o $(S2BIN)/s2_skymask $(S2INC)/s2_skymask.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skydown.o:     $(S2PROG)/s2_skydown.f90 lib
$(S2BIN)/s2_skydown:       $(S2INC)/s2_skydown.o
	$(FC) -o $(S2BIN)/s2_skydown $(S2INC)/s2_skydown.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyup.o:       $(S2PROG)/s2_skyup.f90 lib
$(S2BIN)/s2_skyup:         $(S2INC)/s2_skyup.o
	$(FC) -o $(S2BIN)/s2_skyup $(S2INC)/s2_skyup.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skymultiply.o: $(S2PROG)/s2_skymultiply.f90 lib
$(S2BIN)/s2_skymultiply:   $(S2INC)/s2_skymultiply.o
	$(FC) -o $(S2BIN)/s2_skymultiply $(S2INC)/s2_skymultiply.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyoffset.o: $(S2PROG)/s2_skyoffset.f90 lib
$(S2BIN)/s2_skyoffset:   $(S2INC)/s2_skyoffset.o
	$(FC) -o $(S2BIN)/s2_skyoffset $(S2INC)/s2_skyoffset.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_gcmb.o:        $(S2PROG)/s2_gcmb.f90 lib
$(S2BIN)/s2_gcmb:          $(S2INC)/s2_gcmb.o
	$(FC) -o $(S2BIN)/s2_gcmb $(S2INC)/s2_gcmb.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_grf2.o:        $(S2PROG)/s2_grf2.f90 lib
$(S2BIN)/s2_grf2:          $(S2INC)/s2_grf2.o
	$(FC) -o $(S2BIN)/s2_grf2 $(S2INC)/s2_grf2.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_gcmbcoad.o:        $(S2PROG)/s2_gcmbcoad.f90 lib
$(S2BIN)/s2_gcmbcoad:          $(S2INC)/s2_gcmbcoad.o
	$(FC) -o $(S2BIN)/s2_gcmbcoad $(S2INC)/s2_gcmbcoad.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyrot.o:        $(S2PROG)/s2_skyrot.f90 lib
$(S2BIN)/s2_skyrot:          $(S2INC)/s2_skyrot.o
	$(FC) -o $(S2BIN)/s2_skyrot $(S2INC)/s2_skyrot.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_almrot.o:        $(S2PROG)/s2_almrot.f90 lib
$(S2BIN)/s2_almrot:          $(S2INC)/s2_almrot.o
	$(FC) -o $(S2BIN)/s2_almrot $(S2INC)/s2_almrot.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyorder.o:        $(S2PROG)/s2_skyorder.f90 lib
$(S2BIN)/s2_skyorder:          $(S2INC)/s2_skyorder.o
	$(FC) -o $(S2BIN)/s2_skyorder $(S2INC)/s2_skyorder.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_dlwrite.o:        $(S2PROG)/s2_dlwrite.f90 lib
$(S2BIN)/s2_dlwrite:          $(S2INC)/s2_dlwrite.o
	$(FC) -o $(S2BIN)/s2_dlwrite $(S2INC)/s2_dlwrite.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyerror.o:        $(S2PROG)/s2_skyerror.f90 lib
$(S2BIN)/s2_skyerror:          $(S2INC)/s2_skyerror.o
	$(FC) -o $(S2BIN)/s2_skyerror $(S2INC)/s2_skyerror.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skyrms.o:        $(S2PROG)/s2_skyrms.f90 lib
$(S2BIN)/s2_skyrms:          $(S2INC)/s2_skyrms.o
	$(FC) -o $(S2BIN)/s2_skyrms $(S2INC)/s2_skyrms.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_skythres.o:        $(S2PROG)/s2_skythres.f90 lib
$(S2BIN)/s2_skythres:          $(S2INC)/s2_skythres.o
	$(FC) -o $(S2BIN)/s2_skythres $(S2INC)/s2_skythres.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_ylm.o:        $(S2PROG)/s2_ylm.f90 lib
$(S2BIN)/s2_ylm:          $(S2INC)/s2_ylm.o
	$(FC) -o $(S2BIN)/s2_ylm $(S2INC)/s2_ylm.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_nmask.o:        $(S2PROG)/s2_nmask.f90 lib
$(S2BIN)/s2_nmask:          $(S2INC)/s2_nmask.o
	$(FC) -o $(S2BIN)/s2_nmask $(S2INC)/s2_nmask.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_sky2fsht.o:        $(S2PROG)/s2_sky2fsht.f90 lib
$(S2BIN)/s2_sky2fsht:          $(S2INC)/s2_sky2fsht.o
	$(FC) -o $(S2BIN)/s2_sky2fsht $(S2INC)/s2_sky2fsht.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_map2matmap.o:        $(S2PROG)/s2_map2matmap.f90 lib
$(S2BIN)/s2_map2matmap:          $(S2INC)/s2_map2matmap.o
	$(FC) -o $(S2BIN)/s2_map2matmap $(S2INC)/s2_map2matmap.o \
	$(LDFLAGS) $(PPFLAGS) 

$(S2INC)/s2_alm2matalm.o:        $(S2PROG)/s2_alm2matalm.f90 lib
$(S2BIN)/s2_alm2matalm:          $(S2INC)/s2_alm2matalm.o
	$(FC) -o $(S2BIN)/s2_alm2matalm $(S2INC)/s2_alm2matalm.o \
	$(LDFLAGS) $(PPFLAGS) 
