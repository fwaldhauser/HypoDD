MK	= make
#SUBDIRS	= hypoDD ncsn2pha ph2dt hista2ddsta
SUBDIRS	= hypoDD ncsn2pha ph2dt 

all:
	for d in $(SUBDIRS); do (cd $$d; echo $$d; $(MK) $@); done

clean:
	for d in $(SUBDIRS); do (cd $$d; echo $$d; $(MK) $@); done
