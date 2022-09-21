ciriorbitmax.for, version number can be found at the end of this comment.
c-----------------------------------------------------------------------
c
c Program for computing IRI Fe peak parameters along a satellite orbit
c
c corrections
c-version-mm/dd/yy------------------------------------------------------
c 2000.01 05/07/01 initial version
c 2000.02 07/11/01 line 210: do i=1,100 instead of i=2,100 (K. Tokar)
c 2000.03 28/12/01 output oar(39) for IG12 (R. Conkright, NGDC, NOAA)
c 2000.04 28/10/02 replace TAB/6 blanks, enforce 72/line (D. Simpson)
c 2000.05 02/06/03 Ne(Te) only 300,400; foF1 and hmF1 output corr.
c 2000.06 01/19/05 (.not.jf(20)) instead of (..jf(2)) (G. Schiralli)
c 2005.01 05/06/06 included spread-F (jf(28)) and topside (jf(29)) options
C 2007.00 05/18/07 Release of IRI-2007
c 2007.02 10/31/08 outf(100) -> outf(500), numhei=numstp=500
c 2007.03 02/12/09 added new D-region option (h=-3)
c 2007.11 04/19/10 correct TEC for normal output  [Shunrong Zhang] 
c 2017 to use iri_2016
c
      INTEGER        pad1(6),jdprof(77)
      DIMENSION      outf(20,1000),oar(100),jfi(6)
      DIMENSION      a(20,1000),b(100,1000)
      LOGICAL        jf(50)
      CHARACTER*3    uni(48),sopt
      CHARACTER*4    IMZ(8),MAP,xtex,coorv(2)
      CHARACTER*5    ITEXT(8)
      CHARACTER*6    dopt,pna(48)
      CHARACTER*8    bopt
      CHARACTER*9    topt,pname(7)
      CHARACTER*10   iopt
      CHARACTER*16   f1opt

      CHARACTER*90      outfile1,OUTFILE2,OUTFILE3,OUTFILE4
      CHARACTER*90      outfile5
      CHARACTER*10      cyear,cdoy
      CHARACTER*20     filename
      CHARACTER*20     isr(36),c_st(36) 
      integer*4 time_min
      real*8  tec(37,73),nmf2(37,73),hmf2(37,73), ap_ut(24)

      integer inum_st 
      real  xlon(36),xlat(36)

      data jfi/8,9,13,14,15,16/
      COMMON/const2/icalls,nmono,iyearo,idaynro,rzino,igino,ut0

      icalls=0
      nmono=-1
      iyearo=-1
      idaynro=-1
      rzino=-1
      igino=-1
      ut0=-1

      include './station.fof2.2012.inc'
        call read_ig_rz
        call readapf107
        
        do 6249 i=1,100
6249    oar(i)=-1.0

        do i=1,43
              jf(i)=.true.
        enddo
        jf(2)=.false.           ! no temperatures
        jf(3)=.false.           ! no ion composition
        jf(5)=.false.               ! URSI foF2 model  !! 1_IRI
        jf(6)=.false.               ! Newest ion composition model
        jf(12)=.false.              ! no konsol messages
        jf(21)=.false.           ! ion drift computed
        jf(23)=.false.              ! TTS Te model is standard
c        jf(28)=.true.           ! spread-F computed
        jf(29)=.false.              ! New Topside options
        jf(30)=.false.              ! NeQuick topside
        jf(33)=.false.              ! auroral boundary computed
        jf(35)=.false.              ! no foE storm updating 
        jf(39)=.false.              !  new models for hmF2 (jf(40)=true)
                                    !  Shubin-COSMIC model   


        jm=0
        iut=1
        htec_max=0.0
        hxin=300.0
        jmag=jm
        HEIBEG=hxin
        HEIEND=hxin
        HEISTP=1.      

       i = IARGC()
       if(i.lt.1) then
       write(6,*) 'What Year'
       read(5,*) iyyyy
       write(6,*) 'What Last DOY'
       read(5,*) idoy_end
       else
       call getarg(1,cyear)
       read(cyear,*) iyyyy
       call getarg(2,cdoy)
       read(cdoy,*) idoy_end
       end if

      inum_st =3 

C      iyyyy=2012
C      do idoy=1,366
      do idoy=1,idoy_end
        
       call GETDATE(Iyyyy,IM,ID,idoy)

       ievent=im

       print*,"year/doy/mon/day/hr= "

      write(OUTFILE5,'("IRI."i4.4,i2.2,".in.ap.f107")') iyyyy,ievent

      OPEN(UNIT=50,FILE=OUTFILE5,STATUS='unknown',access='append')
C

       do inhh = 0,23
C       do inhh = 0,0
        print*,iyyyy,idoy ,im,id,inhh
       do inmm = 30,30
          inss = 0     !second
           time_min=inhh*60+inmm
        do ist=1,inum_st     !! # of longitude slices

         along = xlon(ist)
         along = mod(along+720.,360.)
         alati = xlat(ist)

       write(outfile1,'("IRI_0004_2012-",i2.2,"Z_NmF2_",a2,".dat")')
     & ievent,isr(ist)
        write(outfile2,'("IRI_0004_2012-",i2.2,"Z_foF2_",a2,".dat")')
     & ievent,isr(ist)
        write(outfile3,'("IRI_0004_2012-",i2.2,"Z_hmF2_",a2,".dat")')
     & ievent,isr(ist)
        write(outfile4,'("IRI_0004_2012-",i2.2,"Z_TEC_",a2,".dat")')
     & ievent,isr(ist)
C        write(outfile1,'("CTIPE_0001_2013-03Z_NmF2_",a2,".q30.dat")')
C     & isr(ist)
C        write(outfile2,'("CTIPE_0001_2013-03Z_foF2_",a2,".q30.dat")')
C     & isr(ist)


         open(10,file=outfile1,access='append')
         open(20,file=outfile2,access='append')
         open(30,file=outfile3,access='append')
         open(40,file=outfile4,access='append')

C       do ilat=  0,180,5
C          alati = ilat - 90.
C       do ilon = 0,360,5
C          along = ilon

        dhour=inhh*1.0+(inmm+inss/60.)/60.+25.
        mmdd= im*100+id

CCCCCC NmF2 and hmF2
       
        call IRI_SUB(JF,JMAG,ALATI,ALONG,IYYYY,MMDD,DHOUR,
     &    HEIBEG,HEIEND,HEISTP,OUTF,OAR)


CCCCCCC TEC calculation
        htec_max=2000.0   !ccc IRI_0004 2000km with iri-2016
        ivar=1

        dhour=inhh*1.0+(inmm+inss/60.)/60.
        call iri_web(jmag,jf,alati,along,iyyyy,mmdd,1,dhour,
     &          hxin,htec_max,ivar,heibeg,heiend,heistp,a,b)
C     &          hxin,htec_max,ivar,vbeg,vend,vstp,a,b)
CCCCCCccccc

        print*, iyyyy,im,id,inhh,inmm,inss,ALATI,ALONG,
     & oar(2), oar(1), b(37,1)/1E16 
C     & ,tec(ia,io),nmf2(ia,io),hmf2(ia,io)
C        WRITE(30,7118) iyyyy,im,id,inhh,inmm,inss,ALATI,ALONG,
C     & oar(2), oar(1), b(37,1)/1E16 
C     & xlonsta,oar(2), oar(1)/1E6, b(37,1)/1E16
C7118    Format(i4.4,2x,5(i2.2,2x),2F11.2,2X,F10.3,1P2E12.3)

C        enddo !ilon
C        enddo !ilat 
       if(id.eq.1.and.inhh.eq.0.and.inmm.eq.30) then
       write(20,'("# IRI 2016")')
       write(20,'("# ",a30)') c_st(ist)
       write(20,'("# glat: ",f10.2,"N, glon: ",f10.2,"E")')ALATI,ALONG
       write(20,'("#year mon day hh  mm  ss  ms        foF2    ")')
       write(20,'("#                                   [MHz]   ")')
       write(30,'("# IRI 2016")')
       write(30,'("# ",a30)') c_st(ist)
       write(30,'("# glat: ",f10.2,"N, glon: ",f10.2,"E")')ALATI,ALONG
       write(30,'("#year mon day hh  mm  ss  ms        hmF2    ")')
       write(30,'("#                                   [Km]   ")')
       write(30,'("# IRI 2016")')
       write(30,'("# ",a30)') c_st(ist)
       write(40,'("# glat: ",f10.2,"N, glon: ",f10.2,"E")')ALATI,ALONG
       write(40,'("#year mon day hh  mm  ss  ms         TEC    ")')
       write(40,'("#                                   [TECU]   ")')
       endif
       

C       if(xnmf2_2.ne.0.and.xnmf2_1.ne.0.and.imm.eq.30) then
       write(10,'(i4,2x,6(i2.2,2x),1P1E12.4)')
     &       iyyyy,im,id,inhh,inmm,0,0,oar(1)

        write(20,'(i4,2x,6(i2.2,2x),f10.2)')
     &       iyyyy,im,id,inhh,inmm,0,0,sqrt(oar(1)/1.24E10)

        write(30,'(i4,2x,6(i2.2,2x),f10.2)')
     &       iyyyy,im,id,inhh,inmm,0,0,oar(2) 

        write(40,'(i4,2x,6(i2.2,2x),f10.2)')
     &       iyyyy,im,id,inhh,inmm,0,0,b(37,1)/1E16 
C       endif


          close (10)
          close (20)
          close (30)
          close (40)

        enddo !ist (station)
        ap_ut(inhh) = oar(51)
          
C        print*, "upto here", oar(51)

        enddo !inmm
        enddo !inhh
        enddo !idoy

       WRITE(50,102)IYear,IDOY,IM,ID,(ap_ut(ii),ii=0,21,3),
     &  oar(52),oar(41)
        goto 4322
102    FORMAT(i4.4,2x,i3.3,2x,2(i2.2,2x),10F6.0)

4321    print*,'ERROR'
        print*,iyyyy,im,id,inhh,inmm,inss,ALATI,ALONG
        goto 9876
4322    print*,'END'
        print*,iyyyy,im,id,inhh,inmm,inss,ALATI,ALONG

9876           stop
           end

