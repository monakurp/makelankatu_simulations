#!/bin/bash

# ------------------------------------------------------------------------------------------------#
# DEFINE PATHS AND VARIABLES HERE:

# Input filepath:
inpath=/home/monakurp/Helsinki3D

# Script path:
scriptpath=/home/monakurp/makelankatu_simulations/scripts

# booleans:
include_land_use_types=false # include land use types
include_child_lanes=false # include child lane maps

# Horizontal resolutions:
dx_root=9
dx_parent=3
dx_child=1

# Horizontal dimensions:
nx_root=768
ny_root=768
nx_parent=768
ny_parent=768
nx_child=576
ny_child=576

# Origin of the map (N=Northing, E=Easting):
pN=6676070
pE=25497317
pE_root=25496417

# Pivot point (0-1):
rLx=0.5
rLy=0.5

# Wind direction:
wd=270

# Trees: leaf area index (LAI), and the bottom and top height (m) of a reference tree:
lai=12
lai_winter=2.4
tree_bottom=4
tree_top=14

# Fill values
fbyte=-127
fint=-9999
ffloat=-9999.0

# ------------------------------------------------------------------------------------------------#
# Root topography, orography, vegetation and land use:

# Select area
extractDomainFromTile.py -f $inpath/Topography_H3D.npz -fo topo_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
extractDomainFromTile.py -f $inpath/Orography_H3D.npz -fo oro_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
extractDomainFromTile.py -f $inpath/Treeheight_agl_H3D.npz -fo vege_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
extractDomainFromTile.py -f $inpath/Land_cover_types_H3D.npz -fo landuse_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
if $include_land_use_types; then
  extractDomainFromTile.py -f $inpath/Water_type.npz -fo water_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
  extractDomainFromTile.py -f $inpath/Pavements_withMaanpeite.npz -fo pavement_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
  extractDomainFromTile.py -f $inpath/Soil_type.npz -fo soil_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
  extractDomainFromTile.py -f $inpath/Vegetation_type.npz -fo vegetation_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
  extractDomainFromTile.py -f $inpath/Building_type.npz -fo building_type_root -gc -iP $pN $pE_root -dx $dx_root $dx_root -r $rLx $rLy -wd $wd -N $nx_root $ny_root
fi

# Create LAD file
rasterToCanopy3D.py -f vege_root.npz -fo lad_root.npz -m const -l $lai -zr $tree_bottom $tree_top
rasterToCanopy3D.py -f vege_root.npz -fo lad_winter_root.npz -m const -l $lai_winter -zr $tree_bottom $tree_top

# Create building_id map
maskFromRasterTile.py -f landuse_root.npz -fo building_id_root.npz -mv 6


if $include_land_use_types; then

  #replace user-defined values (0) with fill values ($fbyte)
  replaceRasterValues.py -f water_root.npz -fo water_root.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_root $nx_root
  replaceRasterValues.py -f pavement_root.npz -fo pavement_root.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_root $nx_root
  replaceRasterValues.py -f soil_root.npz -fo soil_root.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_root $nx_root
  replaceRasterValues.py -f vegetation_root.npz -fo vegetation_root.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_root $nx_root

  # Approximate undefined soil types as medium porosity (2)
  replaceRasterValues.py -f soil_root.npz -fo soil_root.npz -gt 6 -v 2 -p1 0 0 -p2 $ny_root $nx_root
fi

# Remove orography from topograhy to obtain the building height
$scriptpath/rasterTileSubtraction.py -f topo_root.npz -fs oro_root.npz -fm building_id_root.npz -mvl 6.0 -ntp -fo topo_root.npz

#create building id from topography
replaceRasterValues.py -f topo_root.npz -fo building_id_root.npz -gt 0 -v 1 -p1 0 0 -p2 $ny_root $nx_root

#replace with fill value (-9999) where there are no buildings
replaceRasterValues.py -f topo_root.npz -fo topo_root.npz -lt 0.000000001 -v $ffloat -p1 0 0 -p2 $ny_root $nx_root

#same for building id, but with different fill value
replaceRasterValues.py -f building_id_root.npz -fo building_id_root.npz -lt 1 -v $fint -p1 0 0 -p2 $ny_root $nx_root

# Given a unique building_id value for each building
$scriptpath/fixBuildingID.py -f building_id_root.npz

if $include_land_use_types; then

  #approximate building type as office building built between 1950-2000
  replaceRasterValues.py -f building_id_root.npz -fo building_root.npz -gt 0 -v 5 -p1 0 0 -p2 $ny_root $nx_root

  # Set values samller than 1 to the fill value (-127)
  replaceRasterValues.py -f building_root.npz -fo building_root.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_root $nx_root

  # fill missing surface types according to PIDS
  $scriptpath/fixSurfaceTypes.py -f0 vegetation_root.npz -f1 pavement_root.npz -f2 water_root.npz -f3 building_root.npz -f4 building_type_root.npz -f5 soil_root.npz

  # create surface fractions
  $scriptpath/rastersToSurfaceFraction.py -f0 vegetation_root.npz -f1 pavement_root.npz -f2 water_root.npz -fo surface_fraction_root
fi



# ------------------------------------------------------------------------------------------------#
# Parent topography, orography, vegetation and land use:

# Select area
extractDomainFromTile.py -f $inpath/Topography_H3D.npz -fo topo_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
extractDomainFromTile.py -f $inpath/Orography_H3D.npz -fo oro_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
extractDomainFromTile.py -f $inpath/Treeheight_agl_H3D.npz -fo vege_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
extractDomainFromTile.py -f $inpath/Land_cover_types_H3D.npz -fo landuse_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
if $include_land_use_types; then
  extractDomainFromTile.py -f $inpath/Water_type.npz -fo water_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
  extractDomainFromTile.py -f $inpath/Pavements_withMaanpeite.npz -fo pavement_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
  extractDomainFromTile.py -f $inpath/Soil_type.npz -fo soil_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
  extractDomainFromTile.py -f $inpath/Vegetation_type.npz -fo vegetation_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
  extractDomainFromTile.py -f $inpath/Building_type.npz -fo building_type_parent -gc -iP $pN $pE -dx $dx_parent $dx_parent -r $rLx $rLy -wd $wd -N $nx_parent $ny_parent
fi

# Create LAD file
rasterToCanopy3D.py -f vege_parent.npz -fo lad_parent.npz -m const -l $lai -zr $tree_bottom $tree_top
rasterToCanopy3D.py -f vege_parent.npz -fo lad_winter_parent.npz -m const -l $lai_winter -zr $tree_bottom $tree_top

# Create building_id map
maskFromRasterTile.py -f landuse_parent.npz -fo building_id_parent.npz -mv 6

if $include_land_use_types; then

  #replace user-defined values
  replaceRasterValues.py -f water_parent.npz -fo water_parent.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f pavement_parent.npz -fo pavement_parent.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f soil_parent.npz -fo soil_parent.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f vegetation_parent.npz -fo vegetation_parent.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f soil_parent.npz -fo soil_parent.npz -gt 6 -v 2 -p1 0 0 -p2 $ny_parent $nx_parent

fi

# Remove orography from topograhy to obtain the building height
$scriptpath/rasterTileSubtraction.py -f topo_parent.npz -fs oro_parent.npz -fm building_id_parent.npz -mvl 6.0 -ntp -fo topo_parent.npz

# Given a unique building_id value for each building
$scriptpath/fixBuildingID.py -f building_id_parent.npz

if $include_land_use_types; then

  # same as with root
  replaceRasterValues.py -f topo_parent.npz -fo building_id_parent.npz -gt 0 -v 1 -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f building_id_parent.npz -fo building_parent.npz -gt 0 -v 5 -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f topo_parent.npz -fo topo_parent.npz -lt 0.000000001 -v $ffloat -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f building_id_parent.npz -fo building_id_parent.npz -lt 1 -v $fint -p1 0 0 -p2 $ny_parent $nx_parent
  replaceRasterValues.py -f building_parent.npz -fo building_parent.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_parent $nx_parent

  # fill missing surface types according to PIDS
  $scriptpath/fixSurfaceTypes.py -f0 vegetation_parent.npz -f1 pavement_parent.npz -f2 water_parent.npz -f3 building_parent.npz -f4 building_type_parent.npz -f5 soil_parent.npz

  # create surface fractions
  $scriptpath/rastersToSurfaceFraction.py -f0 vegetation_parent.npz -f1 pavement_parent.npz -f2 water_parent.npz -fo surface_fraction_parent

fi


# ------------------------------------------------------------------------------------------------#
# Child topography, orography, vegetation and land use:

# Select area
extractDomainFromTile.py -f $inpath/Topography_H3D.npz -fo topo_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
extractDomainFromTile.py -f $inpath/Orography_H3D.npz -fo oro_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
extractDomainFromTile.py -f $inpath/Treeheight_agl_H3D.npz -fo vege_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
extractDomainFromTile.py -f $inpath/Land_cover_types_H3D.npz -fo landuse_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
if $include_land_use_types; then
  extractDomainFromTile.py -f $inpath/Water_type.npz -fo water_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
  extractDomainFromTile.py -f $inpath/Pavements_withMaanpeite.npz -fo pavement_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
  extractDomainFromTile.py -f $inpath/Soil_type.npz -fo soil_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
  extractDomainFromTile.py -f $inpath/Vegetation_type.npz -fo vegetation_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
  extractDomainFromTile.py -f $inpath/Building_type.npz -fo building_type_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
fi

# Create LAD file
rasterToCanopy3D.py -f vege_child.npz -fo lad_child.npz -m const -l $lai -zr $tree_bottom $tree_top
rasterToCanopy3D.py -f vege_child.npz -fo lad_winter_child.npz -m const -l $lai_winter -zr $tree_bottom $tree_top

# Create building_id map
maskFromRasterTile.py -f landuse_child.npz -fo building_id_child.npz -mv 6

if $include_land_use_types; then

  #replace user-defined values
  replaceRasterValues.py -f water_child.npz -fo water_child.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f pavement_child.npz -fo pavement_child.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f soil_child.npz -fo soil_child.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f vegetation_child.npz -fo vegetation_child.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f soil_child.npz -fo soil_child.npz -gt 6 -v 2 -p1 0 0 -p2 $ny_child $nx_child

fi

# Remove orography from topograhy to obtain the building height
rasterTileSubtraction.py -f topo_child.npz -fs oro_child.npz -fm building_id_child.npz -mvl 6.0 -ntp -fo topo_child.npz

# Given a unique building_id value for each building
$scriptpath/fixBuildingID.py -f building_id_child.npz

if $include_land_use_types; then

  # same as with root/parent
  replaceRasterValues.py -f topo_child.npz -fo building_id_child.npz -gt 0 -v 1 -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f building_id_child.npz -fo building_child.npz -gt 0 -v 5 -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f topo_child.npz -fo topo_child.npz -lt 0.000000001 -v $ffloat -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f building_id_child.npz -fo building_id_child.npz -lt 1 -v $fint -p1 0 0 -p2 $ny_child $nx_child
  replaceRasterValues.py -f building_child.npz -fo building_child.npz -lt 1 -v $fbyte -p1 0 0 -p2 $ny_child $nx_child

  # fill missing surface types according to PIDS
  $scriptpath/fixSurfaceTypes.py -f0 vegetation_child.npz -f1 pavement_child.npz -f2 water_child.npz -f3 building_child.npz -f4 building_type_child.npz -f5 soil_child.npz

  # create surface fractions
  $scriptpath/rastersToSurfaceFraction.py -f0 vegetation_child.npz -f1 pavement_child.npz -f2 water_child.npz -fo surface_fraction_child

fi

if $include_child_lanes; then
  # With lanes:
  extractDomainFromTile.py -f $inpath/Land_cover_types_H3D.npz -fo lanes_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child

  # With street types:
  extractDomainFromTile.py -f ~/makelankatu_simulations/source_data/roadmap_merged.npz -fo street_types_child -gc -iP $pN $pE -dx $dx_child $dx_child -r $rLx $rLy -wd $wd -N $nx_child $ny_child
fi


# Remove log-files
rm -f .extractDomain*
rm -f .rasterTo*
rm -f .maskFromRaster*
rm -f .replaceRasterValues*
