###################################################
### Scratch space to store cdo terminal scripts ###
###################################################

# run through each model in turn...

######################
### AWI-ESM-1-1-LR ###
######################

# process netcdf files for salinity
cd /Users/home/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Omon/sos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI not bil
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for temperature
cd /Users/home/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Omon/tos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for sea surface height
cd /Users/home/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Omon/zos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for mixed layer depth
cd /Users/home/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/Omon/mlotst/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for sea ice concentration
cd /Users/home/CMIP6/CMIP/AWI/AWI-ESM-1-1-LR/historical/r1i1p1f1/SImon/siconc/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


#################
### INM-CM4-8 ###
#################

# INM data only available as gr1 but regrid for consistency with other data structures
# process netcdf files for temperature
cd /Users/home/CMIP6/CMIP/INM/INM-CM4-8/historical/r1i1p1f1/Omon/tos/gr1/20190530/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for salinity
cd /Users/home/CMIP6/CMIP/INM/INM-CM4-8/historical/r1i1p1f1/Omon/sos/gr1/20190530/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for height
cd /Users/home/CMIP6/CMIP/INM/INM-CM4-8/historical/r1i1p1f1/Omon/zos/gr1/20190530/
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files for sea ice concentration
cd /Users/home/CMIP6/CMIP/INM/INM-CM4-8/historical/r1i1p1f1/SImon/siconc/gr1/20190530/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


##################
### MIROC-ES2L ###
##################

# process netcdf files from MIROC-ES2L for salinity
# use rempabil everywhere else

# temperature
cd /Users/home/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/tos/gn/20190823/
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# salinity
cd /Users/home/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/sos/gn/20190823/
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# height
cd /Users/home/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/Omon/zos/gn/20190823/
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/MIROC/MIROC-ES2L/historical/r1i1p1f2/SImon/siconc/gn/20190823/
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


#####################
### MPI-ESM1-2-LR ###
#####################

# process netcdf files from MPI-ESM1-2-LR for salinity
cd /Users/home/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/sos/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature
cd /Users/home/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/tos/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# height
cd /Users/home/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/zos/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# mixed layer depth
cd /Users/home/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/Omon/mlotst/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/MPI-M/MPI-ESM1-2-LR/historical/r1i1p1f1/SImon/siconc/gn/20190710
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


####################
### IPSL-CM6A-LR ###
####################

# use rmapycon to remove latitudinal line of NAs in middle of Indian ocean
# salinity
cd /Users/home/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/Omon/sos/gn/20180803
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature
cd /Users/home/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/Omon/tos/gn/20180803
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea surface height
cd /Users/home/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/Omon/zos/gn/20180803
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# mixed layer depth
cd /Users/home/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/Omon/mlotst/gn/20180803
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/IPSL/IPSL-CM6A-LR/historical/r1i1p1f1/SImon/siconc/gn/20180803
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

#####################
### ACCESS-ESM1-5 ###
#####################

# salinity
cd /Users/home/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/sos/gn/20191115
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature
cd /Users/home/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/tos/gn/20191115
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea surface height
cd /Users/home/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/zos/gn/20191115
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# mixed layer depth
cd /Users/home/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/Omon/mlotst/gn/20191115
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/SImon/siconc/gn/20200817
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

##################
### MRI-ESM2-0 ###
##################

# mri needs remapycon to avoid weird issues with Antarctic continent
# salinity
cd /Users/home/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/sos/gn/20210311
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature
cd /Users/home/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/tos/gn/20190904
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature sea surface height
cd /Users/home/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/zos/gn/20191205
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# mixed layer depth
cd /Users/home/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/Omon/mlotst/gn/20191205
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/MRI/MRI-ESM2-0/historical/r1i1p1f1/SImon/siconc/gn/20190904
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

##################
### NCAR-CESM2 ###
##################

# salinity
cd /Users/home/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Omon/sos/gn/20190308
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# temperature
cd /Users/home/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Omon/tos/gn/20190308
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea surface height
cd /Users/home/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Omon/zos/gn/20190308
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# mixed layer depth
cd /Users/home/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/Omon/mlotst/gn/20190308
for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# sea ice concentration
cd /Users/home/CMIP6/CMIP/NCAR/CESM2/historical/r1i1p1f1/SImon/siconc/gn/20190308
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# ends
