c hypoDD.inc: Stores parameters that define array dimensions in hypoDD.
c             Modify to fit size of problem and available computer memory.
c	      If 3D raytracing is used, also set model parameters in vel3d.inc.
c Parameter Description:
c MAXEVE:   Max number of events (must be at least the size of the number 
c           of events listed in the event file)
c MAXDATA:  Max number of observations (must be at least the size of the 
c           number of observations).  
c MAXEVE0:  Max number of events used for SVD. If only LSQR is used, 
c           MAXEVE0 can be set to 2 to free up memory. 
c MAXDATA0: Max number of observations used for SVD. If only LSQR is used, 
c           MAXDATA0 can be set to 1 to free up memory. 
c MAXLAY:   Max number of model layers.
c MAXSTA:   Max number of stations.
c MAXCL:    Max number of clusters allowed. 
	integer*4 MAXEVE, MAXLAY, MAXDATA, MAXSTA, MAXEVE0, MAXDATA0
	integer*4 MAXCL

c parameters for large size problems 
      parameter(MAXEVE=   6500,    
     &          MAXDATA=  3000000,
     &          MAXEVE0=  215,	
     &          MAXDATA0= 100000,	
     &          MAXLAY=   50,	
     &          MAXSTA=   400,		
     &          MAXCL=    200)		

cc parameters for small size problems 
c      parameter(MAXEVE=   1000,    
c     &          MAXDATA=  20000,
c     &          MAXEVE0=  20,	
c     &          MAXDATA0= 2000,	
c     &          MAXLAY=   30,	
c     &          MAXSTA=   100,		
c     &          MAXCL=    200)		

