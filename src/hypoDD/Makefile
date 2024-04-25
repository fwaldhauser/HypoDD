CMD	= hypoDD
CC	= gcc
#FC	= g77
FC	= gfortran 
#FC	= f77
SRCS	= $(CMD).f \
	  aprod.f cluster1.f covar.f datum.f \
	  delaz.f delaz2.f direct1.f dist.f dtres.f exist.f \
	  freeunit.f getdata.f getinp.f getinp2.f ifindi.f \
	  indexxi.f juliam.f lsfit_lsqr.f lsfit_svd.f \
	  lsqr.f matmult1.f matmult2.f matmult3.f mdian1.f \
	  normlz.f partials.f ran.f redist.f refract.f \
	  resstat.f scopy.f sdc2.f setorg.f skip.f \
	  snrm2.f sort.f sorti.f sscal.f \
	  svd.f tiddid.f trialsrc.f trimlen.f \
	  ttime.f vmodel.f weighting.f get_vel3d.f \
	  partials_3d.f setorg_3d.f ray_3d.f \
	  partials_1dmm.f partials_1dsr.f nval.f
CSRCS	= datetime_.c hypot_.c 
#CSRCS	= atoangle_.c atoangle.c datetime_.c hypot_.c rpad_.c 
OBJS	= $(SRCS:%.f=%.o) $(CSRCS:%.c=%.o)
INCLDIR	= ../../include
CFLAGS	= -O -I$(INCLDIR)

#LDFLAGS = 
LDFLAGS = -L/usr/lib/gcc/i686-apple-darwin10/4.2.1 -lgcc
listflg = -L -c -O

# Flags for GNU g77 compiler
#FFLAGS	= -O -I$(INCLDIR) -fno-silent -Wall -implicit -g

# Flags for SUN f77 compiler
#FFLAGS	= -O -I$(INCLDIR) -g
FFLAGS	= -O -I$(INCLDIR) 

# Following line needed on HP-UX (hasn't been tested, though).
#LDFLAGS	= +U77

all: $(CMD)

$(CMD): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@

# Extensive lint-like diagnostic listing (SUN f77 only)
hypoDD.lst: $(SRCS)
	f77 -e -Xlist -c $(SRCS)

clean:
	-rm -f $(CMD) *.o core a.out *.fln junk

# Include-file dependencies

dist.o		: $(INCLDIR)/geocoord.inc
redist.o	: $(INCLDIR)/geocoord.inc
sdc2.o		: $(INCLDIR)/geocoord.inc
setorg.o	: $(INCLDIR)/geocoord.inc

cluster1.o	: $(INCLDIR)/hypoDD.inc
dtres.o		: $(INCLDIR)/hypoDD.inc
getdata.o	: $(INCLDIR)/hypoDD.inc
hypoDD.o	: $(INCLDIR)/hypoDD.inc
lsfit_lsqr.o	: $(INCLDIR)/hypoDD.inc
lsfit_svd.o	: $(INCLDIR)/hypoDD.inc
partials.o	: $(INCLDIR)/hypoDD.inc
refract.o	: $(INCLDIR)/hypoDD.inc
resstat.o	: $(INCLDIR)/hypoDD.inc
skip.o		: $(INCLDIR)/hypoDD.inc
tiddid.o	: $(INCLDIR)/hypoDD.inc
trialsrc.o	: $(INCLDIR)/hypoDD.inc
ttime.o		: $(INCLDIR)/hypoDD.inc
vmodel.o	: $(INCLDIR)/hypoDD.inc
weighting.o	: $(INCLDIR)/hypoDD.inc

get_vel3d.o	: $(INCLDIR)/vel3d.inc
partials_3d.o	: $(INCLDIR)/hypoDD.inc $(INCLDIR)/vel3d.inc
setorg_3d.o	: $(INCLDIR)/vel3d.inc
ray_3d.o	: $(INCLDIR)/vel3d.inc

#atoangle_.o	: $(INCLDIR)/compat.h
#atoangle_.o	: $(INCLDIR)/f77types.h
datetime_.o	: $(INCLDIR)/f77types.h
#rpad_.o		: $(INCLDIR)/f77types.h
#sscanf3_.o	: $(INCLDIR)/compat.h
#sscanf3_.o	: $(INCLDIR)/f77types.h
#sscanf4_.o	: $(INCLDIR)/compat.h
#sscanf4_.o	: $(INCLDIR)/f77types.h

