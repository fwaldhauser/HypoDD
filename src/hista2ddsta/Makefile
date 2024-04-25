CMD	= hista2ddsta
#FC	= g77
#FC	= f77
FC	= gfortran 
SRCS	= $(CMD).f
OBJS	= $(SRCS:%.f=%.o)
INCLDIR	= ../../include

# Flags for GNU g77 compiler
#FFLAGS	= -O -I$(INCLDIR) -fno-silent -Wall -implicit

# Flags for SUN f77 compiler
FFLAGS	= -O -I$(INCLDIR)

# Following line needed on HP-UX (hasn't been tested, though).
#LDFLAGS	= +U77

all: $(CMD)

$(CMD): $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $(@F:.o=.f) -o $@

# Extensive lint-like diagnostic listing (SUN f77 only)
hista2ddsta.lst: $(SRCS)
	f77 -e -Xlist -c $(SRCS)

clean:
	-rm -f $(CMD) *.o core a.out *.fln junk
