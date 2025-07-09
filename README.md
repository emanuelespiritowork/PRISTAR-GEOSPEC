# PRISMA-code
PRISMA reading, masking, coregistration, smoothing and regriding

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics
#https://doi.org/10.1016/j.isprsjprs.2024.07.003, https://doi.org/10.5281/zenodo.11547257

#Lorenzo Busatto @ CNR-IREA for prismaread package
#https://github.com/IREA-CNR-MI/prismaread

#Federico Filipponi @ CNR-IGAG for coregistration procedure with Arosics and GDAL

#Lorenzo Parigi @ CNR-IREA for smoothing procedure
# HOW TO USE IT
Download PRISMA_code master github repo into a folder that will be the "PRISMA_code folder". Download the docker image from dockerhub and then enable the docker image using the following command in the Windows terminal:
```cmd
cd /path/
docker load -i /path/NAME_OF_FILE.tar
```
Use the following command in the Windows terminal:
```cmd
docker run --rm -ti -e DISABLE_AUTH=true -p 127.0.0.1:8787:8787 -v C:/your/path/to/PRISMA_code/folder:/space:rw eo/rarosics:latest
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

