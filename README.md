# PRISTAR-GEOSPEC: PRISma Tools for Analysis in R for GEOcoding and SPECtral refinement
<img width="370" height="190" alt="image" src="https://github.com/user-attachments/assets/993f23f3-4940-4789-9ebe-cd34ea30afda" />
<img width="370" height="162" alt="image" src="https://github.com/user-attachments/assets/51eefba9-0e0c-4789-9c15-0aa43b98e723" />

PRISTAR-GEOSPEC tool can help you in:
- importing L0, L1, L2 PRISMA products;
- extracting and cleaning cloud mask for PRISMA;
- creating unique hyperspectral datacube VNIR+SWIR;
- extracting raster of angles;
- extracting central wavelengths and FWHMs;
- extracting and computing sun and sensor angles of PRISMA image (e.g. for your own atmospheric correction procedure); 
- DEM-guided orthoprojection of PRISMA with DEM and Sentinel-2 image using AROSICS and Rational Polynomial Function (RPF);
- coregistration of PRISMA to Sentinel-2 image using AROSICS and gdalwarp;
- spectral smoothing and bad-bands removal;
- regriding and cropping the PRISMA image to a master PRISMA image.

# Credits

#Emanuele Spirito @ CNR-IREA 
#Author

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics and for reading L0 products
#https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busetto @ CNR-IREA for prismaread package 
#https://github.com/IREA-CNR-MI/prismaread
#distributed under GPL-3.0 license

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL, for the maintanance of the Docker Container and any hardware-related solution
#https://github.com/GFZ/arosics
#distributed under Apache-2.0 license
#https://gdal.org/en/stable
#distributed under MIT license

#Lorenzo Parigi @ CNR-IREA for smoothing procedure

#Riccardo Canazza for advice in regrid procedure

#Yulun Wu @ University of Ottawa for PRISMA_angle.py code
#https://github.com/yulunwu8/tmart/blob/main/tmart/AEC/read_PRISMA_vaa.py
#distributed under GPL-3.0 license

# HOW TO INSTALL IT
Download docker from https://www.docker.com/products/docker-desktop/. Install it and install WSL using the Docker procedure. After the installation is finished, restart your PC and open Docker. Download PRISTAR-GEOSPEC master github repo into a folder that will be the _PRISTAR-GEOSPEC folder_. Download the docker image from dockerhub https://hub.docker.com/r/emanuelespiritowork/pristar-geospec/tags. Use the following command in the Windows terminal:
```cmd
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 --memory="24576m" --memory-swap="24576m" -v C:/your/path/to/PRISTAR-GEOSPEC/folder:/space:rw emanuelespiritowork/pristar-geospec:#.#
``` 
where you insert version number where above is #.# and 24576 stands for the amount of RAM to be used (24GB in this case). Then a UNIX console will be opened and you can start your Docker writing:
```cmd
./init
``` 
Open a browser and enter the following URL:
```cmd
localhost:8787
```
An Rstudio Server will be loaded. Go to the right panel and click over the setup:

![image](https://github.com/user-attachments/assets/cce0db0c-e775-450c-8362-9c724885a2c1)

then put inside the box:
```cmd
/space/
```
# HOW TO USE IT
Use case:
1) you want read a L0 PRISMA image. Put _data_SubAcq3_C_SWIR_SURFACE-OBS_Part0_S11.h5_ and _data_SubAcq3_C_VNIR_SURFACE-OBS_Part0_S11.he5_ into _put_PRISMA_he5_and_S2_tif_here_ folder. Then choose a L1 PRISMA image with same or similar view angle to the L0 products and put the .he5 file inside the _put_PRISMA_he5_and_S2_tif_here_ folder. Then open _main.R_ file in Rstudio server and in the procedure_order variable put
```r
   procedure_order <- c("inject","read")
```
2) you want to read a L1 or L2 PRISMA image in .he5 format. Put your PRISMA file in .he5 format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then open _main.R_ file in Rstudio server and in the procedure_order variable put
```r
   procedure_order <- c("read")
```
3) you want to read, generate cloud mask and angle file. Put your PRISMA file in .he5 format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then open _main.R_ file in Rstudio server and in the procedure_order variable put
```r
   procedure_order <- c("read","cloud","atcor")
```
4) you want to read and coregister your PRISMA image to a one-band Sentinel-2 image using gdalwarp. Put your PRISMA file in .he5 format and S2 one-band image in .tif format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then open _main.R_ file in Rstudio server and in the procedure_order variable put
```r
   procedure_order <- c("read","coreg")
```
5) you want to read and orthoproject your PRISMA image. Put your PRISMA file in .he5 format and S2 one-band image in .tif format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then put your DEM image into .tif format into _DEM_ folder. Then open _main.R_ file in Rstudio server and in the procedure_order variable put
```r
   procedure_order <- c("read","ortho")
```
6) you want to read, coregister your PRISMA image and spectral smooth with bad bands removal. Put your PRISMA file in .he5 format and S2 one-band image in .tif format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then open _main.R_ file in Rstudio server and put
```r
full_230_bands <- F
procedure_order <- c("read","coreg","smooth")
```
7) you want to read, coregister your PRISMA image and regrid the PRISMA image to a master image. Put your PRISMA file in .he5 format and S2 one-band image in .tif format in the _put_PRISMA_he5_and_S2_tif_here_ folder. Then put your master PRISMA image in .tif format into _master_image_for_regridding_ folder. Choose type of resampling ("N" for nearest neighbor, "B" for bilinear, "C" for cubic). Then open _main.R_ file in Rstudio server and put
```r
regrid_option <- "N" #that can be either "B" or "C" 
procedure_order <- c("read","coreg","regrid")
```


Click on _source_ and the code will start. 
# WHAT IF A PROBLEM
Ask
