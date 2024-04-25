c--ncsn2pha converts Hypoinverse Y2000 archive phase format into the phase
c  format wanted by the ph2dt program, which is part of the hypoDD process.
c--This program must be run with the input and output files specified
c  on the command line.
c--Both input and output files are organized with an earthquake location
c  header, then a list of phases.
c--Output file can be read either free format or fixed column format.

c--The NCSN convention of naming verical stations is followed when selecting only
c  vertical stations for P phases.

c--The option of converting HVO archive files with HVO station naming (and finding
c  vertical stations for P phases) is invoked with "h" as the third argument:
c  ncsn2pha inputfile outputfile [h]

c 11/2000-, felix@diablo.wr.usgs.gov
c--Modified by Fred Klein 2/2001
c 2015/11/06  changed output format to accomodate negative magnitudes.

      implicit none
      integer*4 k,trimlen,
     & yr,mo,dy,
     & narguments,iargc
      real minlat, minlon, mag, rtime
      real herr, verr, res, dmin
      character line*180, fn0*80 
      character  fn1*30
      logical ex

      character c1*1		!part of station code to test
      character c2*1		!part of station code to test
      character c4*1		!part of station code to test
      integer*4 hvert		!Assumptions for vertical stations 0=NCSN, 1=HVO 
      logical iy2k		!T for y2k format, F for old format

      real pimp			!importance of p arrival
      integer*4 pqual		!assigned weight code of P pick
      character premk*2		!P phase label
      real p_sec		!P arrival time in sec
      real p_time		!P travel time
      real p_wghtr		!P weight derived from assigned code

      real simp			!importance of s arrival
      integer*4 squal		!assigned weight code of s pick
      character sremk*2		!s phase label
      real s_sec		!s arrival time in sec
      real s_time		!s travel time
      real s_wghtr		!s weight derived from assigned code
      character sta*7		!full station code

      integer*4 cuspid, date,hr,min,p_min

      real sec, lon,deglat, deglon, lat
      real depth
      character  str30*30

c-- get input file name:
      narguments = iargc()
      if (narguments.lt.2) stop 'ncsn2pha "inputfile" "outputfile" [h]' 
      call getarg (1,str30)
      inquire (FILE= str30,exist=ex)
      if (.not. ex) stop' >>> ERROR OPENING INPUT DATA FILE.'
      fn0= str30 (1:trimlen(str30))
      open (1,file=fn0, status='unknown') 

c-- open output filename 
      call getarg (2,str30)
      fn1= str30(1:trimlen(str30))
      open (2,file=fn1,status='unknown') 

c-- make NCSN station name assumptions, unless "h" is given as the 3rd argument
c-- "h" means make HVO station assumptions.
      hvert=0
      if (narguments.gt.2) then
        call getarg (3,str30)
        if (str30.eq.'h') hvert=1
      end if

c read header: 
100   read(1,'(a)',end=200)line  		! read header line
      if (line(1:1).eq.'$') goto 100		!skip any shadow lines	

      iy2k= (line(1:2).eq.'19' .or. line(1:2).eq.'20')
      if(.not.iy2k) then
         read (line(1:6),'(i6)') date
         date= date+19000000
         read(line(3:4),'(i2)')mo
         read(line(5:6),'(i2)')dy
         read(line(7:8),'(i2)')hr
         read(line(9:10),'(i2)')min
         read(line(11:14),'(f4.2)')sec
         read(line(15:16),'(f2.0)')deglat
         read(line(18:21),'(f4.2)')minlat
         read(line(22:24),'(f3.0)')deglon
         read(line(26:29),'(f4.2)')minlon
         read(line(30:34),'(f5.2)')depth
         read(line(68:69),'(f2.1)')mag
         read(line(46:49),'(f4.2)')res
         read(line(81:84),'(f4.2)')herr
         read(line(85:88),'(f4.2)')verr
         read(line(129:138),'(i10)')cuspid
      else
         read(line(1:8),'(i8)')date
         read(line(5:6),'(i2)')mo
         read(line(7:8),'(i2)')dy
         read(line(9:10),'(i2)')hr
         read(line(11:12),'(i2)')min
         read(line(13:16),'(f4.2)')sec
         read(line(17:18),'(f2.0)')deglat
         read(line(20:23),'(f4.2)')minlat
         read(line(24:26),'(f3.0)')deglon
         read(line(28:31),'(f4.2)')minlon
         read(line(32:36),'(f5.2)')depth
         read(line(148:150),'(f3.2)')mag
         read(line(49:52),'(f4.2)')res
         read(line(86:89),'(f4.2)')herr
         read(line(90:93),'(f4.2)')verr
         read(line(137:146),'(i10)')cuspid
      endif

      yr= (date/10000)
      lat= deglat+minlat/60
      lon= -deglon-minlon/60
      rtime= hr*1000000 + min*10000 + sec*100

c--Write earthquake location line
      write (2,60) '#',yr,mo,dy,hr,min,sec,lat,lon,depth,mag,
     & herr,verr,res,cuspid

c read data lines and get time
      k= 1
105   read (1,'(a)',end=200) line  
      if (line(1:1).eq.'$') goto 105		!skip any shadow lines	
      if (line(1:3).eq.'   ') goto 100		!terminator line, event end

c--select phase picks: (modified by FWK 2/2001)
      
      if (.not.iy2k) then
        c1=line(4:4)		!4th letter of site code
        c2=line(9:9)		!one letter component code
        c4=line(98:98)		!last letter of 3-letter component code (SEED)
        premk=line(5:6)		!P label
        sremk=line(37:38)	!s label
c--Form station code from 2-letter net & 5-letter site code
        sta=(line(99:100)//line(1:4)//line(95:95))

        read (line(8:8),*,err=105) pqual		!p assigned weight code
        read (line(40:40),*,err=105) squal		!s assigned weight code
        read (line(18:19),'(i2)',err=105) p_min		!min for both p & s
        read (line(20:24),'(f5.2)',err=105) p_sec	!P arrival time
        read (line(32:36),'(f5.2)',err=105) s_sec	!s arrival time
        read (line(83:86),'(f4.3)',err=105) pimp 	!importance of P arr
        read (line(87:90),'(f4.3)',err=105) simp 	!importance of s arr

      else					!!Y2K
        c1=line(4:4)		!4th letter of site code
        c2=line(9:9)		!one letter component code
        c4=line(12:12)		!last letter of 3-letter component code (SEED)
        premk=line(14:15)	!P label
        sremk=line(47:48)	!s label
c--Form station code from 2-letter net & 5-letter site code
        sta=(line(6:7)//line(1:5))

c--Do not accept any hand timed data before 1954 (NCSN only)
        if (line(109:109).eq.'H' .and. yr.lt.1954) goto 105 

        read (line(17:17),*,err=105) pqual		!p assigned weight code
        read (line(50:50),*,err=105) squal		!s assigned weight code
        read (line(28:29),'(i2)',err=105) p_min		!min for both p & s
        read (line(30:34),'(f5.2)',err=105) p_sec	!P arrival time
        read (line(42:46),'(f5.2)',err=105) s_sec	!s arrival time
        read (line(101:104),'(f4.3)',err=105) pimp 	!importance of P arr
        read (line(105:108),'(f4.3)',err=105) simp 	!importance of s arr
      endif

c----------------------------------------------------
c--Process a P time if it is valid and weighted
      if (pqual.lt.4 .and. premk.ne.'  ') then

c--P must be from a vertical component. 
c--NCSN has a Z in 3rd letter of 3-letter comp code
        if ((hvert.eq.0 .and. c4.eq.'Z') .or.
c--HVO has a V in 4th letter of site code
     2    (hvert.eq.1 .and. c1.eq.'V')) then

c--Correct weights transformation   (NCSN standard weights)    FWK
          if(pqual.eq.0) then
             p_wghtr= 1.0 
          elseif(pqual.eq.1) then
             p_wghtr= 0.5 
          elseif(pqual.eq.2) then
             p_wghtr= 0.2 
          elseif(pqual.eq.3) then
             p_wghtr= 0.1 
          else
             p_wghtr= 0.0
          endif

c if importance is larger than 0.5:
          if (pimp.ge.0.5) p_wghtr= -p_wghtr

c get travel time
          if (p_min.lt.min) then
            dmin= (p_min+60) - min
          else
            dmin= p_min - min
          endif
          p_time= dmin*60. + (p_sec-sec) 
          write (2,61) sta, p_time, p_wghtr, 'P'
        end if
      end if

c----------------------------------------------------
c--Process an s time if it is valid and weighted
      if (squal.lt.4 .and. sremk.ne.'  ') then

c--S must be from a horizontal component. 
c--NCSN has a Z in 3rd letter of 3-letter comp code
        if ((hvert.eq.0 .and. c4.ne.'Z') .or.
c--HVO has a V in 4th letter of site code
     2    (hvert.eq.1 .and. c1.ne.'V')) then

c--Correct weights transformation   (NCSN standard weights)    FWK
          if(squal.eq.0) then
             s_wghtr= 1.0 
          elseif(squal.eq.1) then
             s_wghtr= 0.5 
          elseif(squal.eq.2) then
             s_wghtr= 0.2 
          elseif(squal.eq.3) then
             s_wghtr= 0.1 
          else
             s_wghtr= 0.0
          endif

c if importance is larger than 0.5:
          if (simp.ge.0.5) s_wghtr= -s_wghtr

c get travel time
          if (p_min.lt.min) then
            dmin= (p_min+60) - min
          else
            dmin= p_min - min
          endif
          s_time= dmin*60. + (s_sec-sec) 
          write (2,61) sta, s_time, s_wghtr, 'S'
        end if
      end if

c----------------------------------------------------
      goto 105

c--Earthquake location format. Expanded for 10-digit ID numbers
c  and large but occasional error numbers.  FWK
60    format (a1,i5,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2,1x,
     & f8.4,1x,f9.4,1x,
     & f7.2,f6.2,f6.2,f6.2,f6.2,1x,i10)
cfw151106     & f7.2,f5.2,f6.2,f6.2,f6.2,1x,i10)

c--Phase format
61    format (a7,f10.3,f8.3,3x,a1)
      
200   close(1)
      close(2)
      end



      integer function TRIMLEN(t)
c------------------------------------------------------------------
c     Author:  Urs Kradolfer, June 1986
c     Call:    nc=TRIMLEN(char)
c
c          --> nc says, how many characters the input-string has
c              (ignoring trailing blanks!).
c
      implicit none
c
      character t*(*)
      do 1 trimlen=LEN(t),1,-1
    1    if(t(trimlen:trimlen).ne.' ')RETURN
      trimlen=1
      end ! of integer function trimlen
c
