c ph2dt.inc: Stores parameters that define array dimensions in ph2dt.
c            Modify to fit size of problem and available computer memory.
c Parameter Description:
c MEV:   Max number of events.
c MSTA:  Max number of stations.
c MOBS:  Max number of phases (P&S) per eventer event.

	integer	MEV, MSTA, MOBS

	parameter(MEV=    12000,    
     &		MSTA=   2600,		
     &		MOBS=   1000)

