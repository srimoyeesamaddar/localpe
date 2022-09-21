#f77 -o iri_tec_rgnl  iri_tec_2.5x5x5.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f
#f77 -o iri_3D iriNeTeTiTn_3D.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f
#f77 -o iri_2D iri_2D.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f

pgf95 -o  iri_2D_2013_mar   iri_2D.2013.mar.f     irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  
pgf95 -o  iri_2D_1yr_2012   iri_2D.1yr.2012.f     irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  

#pgf95 -o  iri_2D_user.exe  iri_2D.user.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.2D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib

#pgf95 -o  iri_2D.exe   iri_2D.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.2D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib

#pgf95 -o  iri_3D.exe   iri_3D.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.2D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib
#pgf95 -o iri iritest.for irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

#pgf95 -o  test.iri_2D_user.exe  tmp.iri_2D.f     irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.2D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib


#pgf95 -o  iri_2D_user.exe  iri_2D.user.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.2D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib


#pgf95 -o  iri_3D.exe   iri_3D.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f  netcdf-write.3D.f netcdf.f -lnetcdf -lm -I/public/Linux_64/netcdf-4.1.2_no_underscore/include -L/public/Linux_64/netcdf-4.1.2_no_underscore/lib
#f77 -o   iri_2D_fof2_at_1_location   iri_2D_fof2_at_1_location.f    irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f

#f77 -o   iri_2D_RoR_at_1_location   iri_2D_RoR_at_1_location.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f

#f77 -o  anna_iri_2D anna.iri_nmhmF2TEC.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for dateJVE.f
#f77 -o  iri_gps_tec_D iri_gps_tec_D.f  irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
#g77 -o  iri_gps_tec_D.1yr  iri_gps_tec_D.1yr.f  irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o   iri_gps_tec_global   iri_gps_tec_global.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iri_nmhmF2_global iri_nmhmF2_global.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
#f77 -o iri_nmhmF2_global.1yr iri_nmhmF2_global.1yr.f   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iri.temp  iritemp.for   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iri.tec.gps  iri_tec.gps.for   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iri2.eden.champ  iri.eden.champ.for   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
#f77 -o iri2.NeTeTi300  iriNeTeTi300_db.for   irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
#f77 -o iriorbitmax iriorbitmax.for  irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iriorbit.test iriorbit.for  irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
#f77 -o iriVdrift.test iriVdrift.for  irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for
