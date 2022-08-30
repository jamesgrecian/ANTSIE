###################################################
### Scratch space to store cdo terminal scripts ###
###################################################

################
### Salinity ###
################

# process netcdf files from AWI-ESM-1-1-LR for salinity
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/AWI/AWI-ESM-1-1-LR/lgm/r1i1p1f1/Omon/sos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MIROC-ES2L for salinity
# use rempabil everywhere else
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MIROC/MIROC-ES2L/lgm/r1i1p1f2/Omon/sos/gn/20191002/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MPI-ESM1-2-LR for salinity
cd ~/Dropbox/git_projects/ANTSIE/data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/sos/gn/20190710/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


###################
### Temperature ###
###################

# process netcdf files from AWI-ESM-1-1-LR for temperature
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/AWI/AWI-ESM-1-1-LR/lgm/r1i1p1f1/Omon/tos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from INM-CM4-8 for temperature
# use rempabil everywhere else
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/INM/INM-CM4-8/lgm/r1i1p1f1/Omon/tos/gr1/20201112/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MIROC-ES2L for temperature
# use rempabil everywhere else
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MIROC/MIROC-ES2L/lgm/r1i1p1f2/Omon/tos/gn/20191002/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


# process netcdf files from MPI-ESM1-2-LR for temperature
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/tos/gn/20190710/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

##########################
### Sea Surface Height ###
##########################

# process netcdf files from AWI-ESM-1-1-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/AWI/AWI-ESM-1-1-LR/lgm/r1i1p1f1/Omon/zos/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MIROC-ES2L for sea surface height
# use rempabil everywhere else
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MIROC/MIROC-ES2L/lgm/r1i1p1f2/Omon/zos/gn/20191002/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


# process netcdf files from MPI-ESM1-2-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/zos/gn/20190710/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

#############################
### Sea Ice Concentration ###
#############################

# process netcdf files from AWI-ESM-1-1-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/AWI/AWI-ESM-1-1-LR/lgm/r1i1p1f1/SImon/siconc/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MIROC-ES2L for sea surface height
# use rempabil everywhere else
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MIROC/MIROC-ES2L/lgm/r1i1p1f2/SImon/siconc/gn/20191002/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done


# process netcdf files from MPI-ESM1-2-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/SImon/siconc/gn/20190710/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

#########################
### Mixed Layer Depth ###
#########################

# process netcdf files from AWI-ESM-1-1-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/AWI/AWI-ESM-1-1-LR/lgm/r1i1p1f1/Omon/mlotst/gn/20200212/

# loop through .nc files
# note use remapycon for AWI
for i in *.nc;
do
echo $i
cdo remapycon,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done

# process netcdf files from MPI-ESM1-2-LR for sea surface height
cd /Volumes/Nifty\ 128/ANTSIE\ data/CMIP6/PMIP/MPI-M/MPI-ESM1-2-LR/lgm/r1i1p1f1/Omon/mlotst/gn/20190710/

for i in *.nc;
do
echo $i
cdo remapbil,global_1  "$i" "${i%.nc}_bil_1x1.nc"
done



