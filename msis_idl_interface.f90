program msis_idl_interface  !Called by IDL

!this code is the fortran wrapper that calls the input idl files and
!calls the gtd8 subroutine and
! stores the output to a file
!NOTE: Convert this into a sharable object as below:
!Intel (Windows) and gfortran (Linux). The validated gfortran command is:


!"C:\cygwin\bin\gfortran" -ffree-line-length-200 - msis_idl_interface.f90 hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 msis1.97_test.f90 -o msis_idl_interface.exe

!______________________________________________________________________________________________________________________________________________

!Inputs
!hardwiring number of altitudes for now #alt=400

       integer nalt,i
       real*4 iyd, utsec, lat, lon, stl, f107a, f107, ap, mass,alt(400)


! Outputs-initailize the input varibles for the output file
       
       real*4 d(9),t(2)
       real*4 oden(400) ,o2den(400) ,n2den(400), temp(400)
!______________________________________________________________________________________________________________________________________________

!Open input and outfiles for read and write operations

       open(1, file='/home/srimoyee/Desktop/nrl_files/sav_files/msis_inputs.dat')  !input file
       open(2, file='/home/srimoyee/Desktop/nrl_files/sav_files/msis_outputs.dat')  ! output file

!______________________________________________________________________________________________________________________________________________
      
!Read the input file
       read(1,*)        !skipping the header
       read(1,*) iyd
       read(1,*) utsec
       read(1,*) nalt
      
       read(1,*)  alt(1:nalt)
       read(1,*)  lat
       read(1,*)  lon
       read(1,*) stl
       read(1,*) f107a
       read(1,*) f107
       read(1,*) ap
       read(1,*) mass
       
   !    print *, iyd,utsec,nalt,  lat, lon, stl, f107a, f107, ap, mass
   !    print *,alt
      
!______________________________________________________________________________________________________________________________________________

! write the header to the output file
       write(2,*) 'O(cm-3)', 'O2 (cm-3)' ,'N2(cm-3)','Temperature (K)'
  
  !call the main msis subroutine
       do i=1,nalt
          CALL gtd8(iyd,utsec,alt(i),lat,lon,stl,f107a,f107,ap,mass,d,t)
      
!       MSIS OUTPUT VARIABLES:
!        D(1) - HE NUMBER DENSITY(CM-3)
!        D(2) - O NUMBER DENSITY(CM-3)
!        D(3) - N2 NUMBER DENSITY(CM-3)
!        D(4) - O2 NUMBER DENSITY(CM-3)
!        D(5) - AR NUMBER DENSITY(CM-3)
!        D(6) - TOTAL MASS DENSITY (including anomalous oxygen)(GM/CM3)
!        D(7) - H NUMBER DENSITY(CM-3)
!        D(8) - N NUMBER DENSITY(CM-3)
!        D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
!        T(1) - EXOSPHERIC TEMPERATURE(Not available in 1.97)
!        T(2) - TEMPERATURE AT ALTITUDE(K)
      
             oden(i)= d(2)
             o2den(i)= d(4)
             n2den(i)= d(3)
             temp(i)= t(2)
        
        
!  Write the outputs in a text file to be accessed by idl
             write(2,*) oden(i), o2den(i), n2den(i), temp(i)
             print *, alt(i),temp(i)
      
          end do
   
!______________________________________________________________________________________________________________________________________________


!Close the files
   
      close(1)
      close(2)
                       
                                                                 

       END !subroutine msis_idl_interface
 
 
 
 
