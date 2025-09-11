# PRISMA-code
PRISMA reading, masking, coregistration, smoothing and regriding

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics
#https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busetto @ CNR-IREA for prismaread package 
#https://github.com/IREA-CNR-MI/prismaread
#distributed under GPL-3.0 license

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL
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
Download docker from https://www.docker.com/products/docker-desktop/. Download PRISMA_code master github repo into a folder that will be the _PRISMA_code folder_. Download the docker image from dockerhub https://hub.docker.com/r/emanuelespiritowork/pristar-geospec/tags. Use the following command in the Windows terminal:
```cmd
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 -v C:/your/path/to/PRISMA_code/folder:/space:rw emanuelespiritowork/pristar-geospec:#.#
``` 
where you insert version number where above is #.#. Open a browser and enter the following URL:
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
Open the _main.R_ file in the Rstudio server and put the PRISMA file in .he5 format and the S2 reference file in the .tif format in the folder _put_PRISMA_he5_and_S2_tif_here_. Click on _source_ and the code will start. 
# WHAT IF A PROBLEM
Ask
