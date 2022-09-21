===============================================================================
MSIS® (NRL-SOF-014-1) SOFTWARE

MSIS® is a registered trademark of the Government of the United States of 
America, as represented by the Secretary of the Navy. Unauthorized use of 
the trademark is prohibited. 

The MSIS® Software (hereinafter Software) is property of the United States 
Government, as represented by the Secretary of the Navy. The Government of 
the United States of America, as represented by the Secretary of the Navy, 
herein grants a non-exclusive, non-transferable license to the Software for 
academic, non-commercial, purposes only. A user of the Software shall not: 
(i) use the Software for any non-academic, commercial purposes, (ii) make 
any modification or improvement to the Software, (iii) disseminate the 
Software or any supporting data to any other person or entity who will use 
the Software for any non-academic, commercial purposes, or (iv) copy the 
Software or any documentation related thereto except for (a) distribution 
among the user’s personal computer systems, archival, or emergency repair 
purposes, or (b) distribution for non-commercial, academic purposes, without 
first obtaining the written consent of IP Counsel for the Naval Research 
Laboratory. 

As the owner of MSIS®, the United States, the United States Department of 
Defense, and their employees: (1) Disclaim any warranties, express, or 
implied, including but not limited to any implied warranties of 
merchantability, fitness for a particular purpose, title or non-infringement, 
(2) Do not assume any legal liability or responsibility for the accuracy, 
completeness, or usefulness of the software, (3) Do not represent that use of 
the software would not infringe privately owned rights, (4) Do not warrant 

that the software will function uninterrupted, that is error-free or that any 
errors will be corrected.

BY USING THIS SOFTWARE YOU ARE AGREEING TO THE ABOVE TERMS AND CONDITIONS.  
===============================================================================

NRLMSIS 1.97 Neutral Atmosphere Empirical Model

This software package is a beta version of an upgrade to the
NRLMSISE-00 model of atmospheric temperature and species densities.
Please do not redistribute without contacting the authors.

VERSION HISTORY
  08 MAR 19 Version 1.97 (Beta version)

AUTHORS
  Douglas Drob (douglas.drob@nrl.navy.mil)
  John Emmert (john.emmert@nrl.navy.mil)
  
PACKAGE CONTENTS
  readme.txt            This file
  msis1.97_test.f90     Test driver
  msis1.97_gtd8.f90     Source subroutine to evaluate the model using
                          the legacy interface
  pymsis.f90            Subroutines to initialize the model, set
                          switches, and evaluate the model using a new
                          interface (callable from python via f2py.py)
  constants.f90         Module containing model constants
  hgh2gph.f90           Subroutines to convert between geometric height
                          and geopotential height
  msis.f90              Subroutines to calculate horizontal expansion
                          functions
  tfn.f90               Subroutines to calculate vertical temperature
                          profile
  dfn.f90               Subroutines to calculate vertical density
                          profile
  msis1.97.bin          Binary data file containing model parameters
  msis1.97_test_in.txt  ASCII file containing input for test driver.
  msis1.97_test_out0.txt  ASCII file containing expected output of test
                            driver.
 
COMPILER NOTES
  The model package was tested using the following Fortran compilers:
  Intel (Windows) and gfortran (Linux). The validated gfortran command
  is:
    gfortran -ffree-line-length-200 -o test.exe hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 msis1.97_test.f90  

 "C:\cygwin64\etc\setup\cygwin32-gcc-fortran" -ffree-line-length-200 -o test.exe hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 msis1.97_test.f90  

    "C:\cygwin\bin\gfortran" -ffree-line-length-200 -o msis_idl_interface.exe msis_idl_interface.f90 hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 msis1.97_test.f90  
"


    "C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin\i686-w64-mingw32-gfortran" -ffree-line-length-200 -o msis_idl_interface.exe msis_idl_interface.f90 hgt2gph.f90 constants.f90 msis.f90 tfn.f90 dfn.f90 pymsis.f90 msis1.97_gtd8.f90 
"

RELEASE NOTES: MODEL FORMULATION
  Major changes to the NRLMSISE-00 formulation include:
  - The transition from a fully mixed atmosphere to diffusive separation
    is now represented via height-dependent effective mass for each
    species.
  - The temperature profile is now C2 continuous.
  - A global geopotential height function is now used internally.
  - Atomic oxygen now extends down to 50 km; below 85 km, the O profile
    is represented by cubic B splines decoupled from temperature.
  - Thermal diffusion is no longer applied to any species.
 
 RELEASE NOTES: PARAMETER ESTIMATION
  The parameters of the model were tuned to extensive new data in the
  troposphere, stratosphere, and mesosphere, and to NRLMSISE-00 in
  the thermosphere (i.e., new lower and middle atmosphere, old
  thermosphere). Additional tuning to new thermospheric data is planned
  for the v2.0 release.
