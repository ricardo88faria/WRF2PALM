# WRF2PALM
This is a set of PALM tools to create topography static driver netcdf, land use that contains all static information and large-scale forcing data, provided from WRF, as input for PALM following PALM Input Data Standard (PIDS) v1.9.
There is also a simple check PALM divisible number of processors for the number of grids script.

## Usage:


1. Topography (create_static_topo):
To use create_static_topo.py, user need to grab a .asc file containing topographic information inside raw_topo folder. User must setup, dx, dy, dz, nz, north, south, east and west, PALM configurations and area limits inside create_static_topo.py. After script runs a PALM compatible static file is created containing topography. A .cgf file, that contains PALM domain configuration (dx, dy, dz, nx, ny, nz, north, south, east and west), is also created to use with the remain scripts

2. Land Use (create_static_landuse):
create_static_landuse.py script load a .tif file inside raw_landuse, use .cgf and static file created by create_static_topo.py. After run this script, static file will be modified with land use data converted to PALM compatible data.

3. WRF forcing (create_static_topo):
In order to create WRF forcing create_dynamic.py script is run. WRF outputs must be placed inside raw_forcing folder. Script load north, south, east and west inside .cgf file, extract boundary conditions in the area set by user and interpolate to PALM configurations dx, dy, dz, nx, ny, nz created by create_static_topo.py script.
If user didnt run create_static_topo script, a .cgf file must be setup with dx, dy, dz, nx, ny, nz, north, south, east and west fields.


**Contacts:**

<ricardo88faria@gmail.com>
