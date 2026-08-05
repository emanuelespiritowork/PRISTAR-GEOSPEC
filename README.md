# PRISTAR-GEOSPEC: PRISma Tools for Analysis in R for GEOcoding and SPECtral refinement
<img width="1808" height="214" alt="logo_CNR_esteso" src="https://github.com/user-attachments/assets/e9fcb8e5-e10f-47df-82ce-ada344309316" >
<img width="200" height="80" alt="image" src="https://github.com/user-attachments/assets/51eefba9-0e0c-4789-9c15-0aa43b98e723" >
<img width="150" height="80" alt="image" src="https://github.com/user-attachments/assets/293419e9-1c48-4afc-868a-9fb82d6e3f33" >
<img width="220" height="80" alt="image" src="https://github.com/user-attachments/assets/ba5df6d6-2df0-4114-9c9f-2258bff6ab6f" />
<img width="230" height="90" alt="image" src="https://github.com/user-attachments/assets/4811e70e-1a9e-4cba-b5d3-cf6339ed994a" />


---

PRISMA (PRecursore IperSpettrale della Missione Applicativa) is a pilot hyperspectral satellite launched by ASI (Agenzia Spaziale Italiana) in 2019. Its main research pillars are agriculture, forest, waters, climate change, and raw material exploration. In the context of PRISMA data processing, three issues arise when dealing with L1 (Top of Atmosphere Radiance) and L2 (Bottom of Atmosphere Reflectance) products. First, the products are delivered in HDF format, which is not commonly ready-to-use in a geospatial software, and hyperspectral data is stored into two different stacks, Visibile to Near InfraRed and ShortWave InfraRed, to be merged for further processing. Secondly, PRISMA geocoding accuracy is not within the pixel size (30 m), and then requires a refinement. Thirdly, users may need to perform an ad hoc atmospheric correction for particular targets (e.g. water) and L1 is not geocoded. Available PRISMA tools, importing its format are ENVI-toolkit (NV5 Geospatial Software, Inc.), PRISMA-toolbox (Planetek) and EnMap-Box (GFZ), but only the former – that is commercial - provides for geocoding refinement and no geocoding for L1. To fill this lack with open SW, we realized a seamless procedure starting from standard HDF files to build the L1 geocoded and the L2 geocoding-refined products, in addition, regridding and smoothing routines are offered. The tool is based on open source libraries: prismaread for data conversion, cube merging and basic metadata gathering, and AROSICS and gdal for the geocoding (refinement of L2 or process of L1). For L1 information useful for atmospheric correction (e.g. sun and sensors view angles, band centers, etc.) are also extracted from HDF metadata. The regridding step is added both to L1 and L2 image to enable for multi-temporal image stacks. Last step is the spectral smoothing that is straightforwardly applied to the L2 process chain while it is applied to L1 chain only after atmospheric correction. Since the workflow takes advantages of different libraries and languages (R, Python) an Rstudio Server docker has been created together with a GitHub code repository, to create an easy-to-use distribution which do not require to solve all the dependencies of all libraries used.

---


PRISTAR-GEOSPEC tool can help you in:
- importing L0, L1, L2 PRISMA products downloaded from ASI portal;
- extracting and cleaning cloud mask for PRISMA;
- creating unique hyperspectral datacube VNIR+SWIR;
- extracting raster of angles;
- extracting central wavelengths and FWHMs;
- extracting and computing sun and sensor angles of PRISMA image (e.g. for your own atmospheric correction procedure);
- using a built-in atmospheric corrector (JPL-NASA ISOFIT) to process your PRISMA L1 products;
- DEM-guided orthoprojection of PRISMA with DEM and Sentinel-2 image using AROSICS and Rational Polynomial Function (RPF);
- coregistration of PRISMA to Sentinel-2 image using AROSICS and gdalwarp. For images where GNSS has gone wrong, this procedure is helped with a rigid shift that you can do with a built-in tool;
- spectral smoothing and bad-bands removal;
- regriding and cropping the PRISMA image to a master PRISMA image so as you can create time series stack of PRISMA images;
- keeping all the useful metadata (central wavelengths, FWHMs, band names, cloud mask) of the image;
- having each intermediate step of all the processing done over the original ASI PRISMA "standard" image in a standardized folder structure and with a coherent naming convention;
- apply all this workflow to an entire archive of time series of PRISMA images over the same area without much effort in organising the data thanks to the already-prepared folder structure.

# Credits

#Emanuele Spirito @ CNR-IREA (first author)

#Lorenzo Parigi @ CNR-IREA (second author) for writing the smoothing procedure, helping with three dockers composition, maintenance and docker networks, finding the ISOFIT integration best practice

#Federico Filipponi @ CNR-IGAG for his coregistration procedure made with Arosics and GDAL, for the maintenance of the first version of the Docker Container and any hardware-related solution

#Giandomenico De Luca @ CNR-IBE for advice on versions of GDAL and Arosics and for injection of L0 products

#Riccardo Canazza @ CNR-IREA for advice in regrid procedure

#Gabriele Candiani @ CNR-IREA for advisory in naming convention, metadata quality and visual check of results

#Rodolfo Ceriani @ UNIMI-DISAA for user-advisory

#Mirco Boschetti @ CNR-IREA for accelerating the path

#Monica Pepe @ CNR-IREA for guiding the whole procedure, scientific knowledge and coordination

# Other references

#AROSICS: https://github.com/GFZ/arosics distributed under Apache-2.0 license

#GDAL: https://gdal.org/en/stable distributed under MIT license

#ISOFIT: https://isofit.github.io/isofit/4.1.0/ distributed under Apache-2.0 license

#Lorenzo Busetto @ CNR-IREA for prismaread package https://github.com/IREA-CNR-MI/prismaread distributed under GPL-3.0 license

#Giandomenico De Luca @ CNR-IBE https://doi.org/10.1016/j.isprsjprs.2024.07.003
https://doi.org/10.5281/zenodo.11547257

#Yulun Wu @ Ottawa University https://github.com/yulunwu8/tmart/blob/main/tmart/AEC/read_PRISMA_vaa.py distributed under GPL-3.0 license

#Spirito et al., 2026, Trends in Earth Observation, Volume IV, https://aitonline.org/wp-content/uploads/2026/07/Smart-Earth-Observation-for-a-Sustainable-Future.pdf distributed under Creative Commons AttributionNo Derivatives 4.0 International License

# HOW TO INSTALL IT
Download docker from https://www.docker.com/products/docker-desktop/. Install it and install WSL using the Docker procedure. After the installation is finished, restart your PC and open Docker. Download PRISTAR-GEOSPEC Release v0.9.0-beta and after unzipping this will be your `config_folder`. The maximum RAM used should be increased in the Docker configuration file. Create a file in `C:\Users\<your_users>\.wslconfig` and write inside:
```txt
[wsl2]
memory=24GB
```
or any amount of memory you should maximally give to any Docker. Then in the terminal restart the wsl:
```cmd
wsl --shutdown
```
The PRISTAR_GEOSPEC configuration file for the Docker is the `config_folder/.env` file. Set the amount of RAM you want to limit each docker. Then open the terminal inside that folder (open a terminal and use the `cd` command to put the directory of the `config_folder`) and run:
```cmd
docker compose up -d
``` 
This will check your docker containers and if you are running it for the first time it will download all the dockers needed (namely AROSICS docker, ISOFIT docker and RStudio docker). Then open a browser (Chrome, Edge, Brave, Firefox, ...) and enter the following URL:
```cmd
localhost:8787
```
An Rstudio Server will be loaded. Go to the right panel and click over the setup as in the screenshot:

![image](https://github.com/user-attachments/assets/cce0db0c-e775-450c-8362-9c724885a2c1)

then put inside the box:
```cmd
/config_folder/
```
# WHERE TO PUT YOUR DATA
- (optional) in the `config_folder/DEM` folder you should put the DTM, the aspect and the slope rasters if you want to use the `ortho` feature and the `isofit` atmospheric correction. Each of these files should contain in the name respectively `dtm`, `aspect`, `slope`. In the picture you see the expected `config_folder/DEM` folder content:

<img width="628" height="119" alt="image" src="https://github.com/user-attachments/assets/887388db-3b7b-425a-a34a-62c264f429c6" />

- (optional) in the `config_folder/master_image_for_regridding` folder you should put the master image if you want to use the `regrid` feature and the `crop` feature. The image should be in `.tif` format.  In the picture you see the expected `config_folder/master_image_for_regridding` folder content:

<img width="787" height="75" alt="image" src="https://github.com/user-attachments/assets/d6207eba-d865-4ce6-b8ac-e759f4e1bf78" />

- (mandatory) in the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder you should put a folder for each image you want to process. Inside each of these folders there should be the PRISMA `.he5` file with the original name from ASI portal (starting with `PRS`). Then if you also want to use `coreg` or `ortho` feature then you also need to put here a Sentinel-2 single-band image with the band you want to use to do the coregistration (see `PRS_band_for_coreg` to match to the PRISMA image band used for coregistration). The S2 image should have in its name `S2` or `s2` and should be a `.tif` image. In the picture you see the expected `config_folder/put_PRISMA_he5_and_S2_tif_here` folder content:

<img width="679" height="141" alt="image" src="https://github.com/user-attachments/assets/bc4736b6-c45d-46ae-ba9a-64fa2d575b22" />

An example of content of the `config_folder/put_PRISMA_he5_and_S2_tif_here/PRS_L1_STD_OFFL_20240928` for PRISMA image reading and coregistration should be:

<img width="867" height="55" alt="image" src="https://github.com/user-attachments/assets/ef8ebe49-c1b5-463c-a2c4-763d62213240" />



# HOW TO USE IT
You will see a list of files. Click on `main.R`. Here you will see the script configuration file where you can manage the processing of the images. Here is the list of the parameters with their explaination:

- `regrid_option`: when using the `regrid` step this will be the regriding method, namely `N` for nearest neighbour, `C` for cubic-spline and `B` for bilinear. If you don't know the default `N` will be a conservative choice;
- `full_230_bands`: when using the `smooth` step this will interpolate the bad bands (if `TRUE`) or will just remove the bad bands (if `FALSE`);
- `PRS_band_for_coreg`: when using the `coreg` step this will be the band number to use for coregistration procedure;
- `shift`: when using the `coreg` step this will enable the rigid shift of the image (if `TRUE`);
- `shift_x`: when `shift` is `TRUE`, this will be the meters of the shift along x-axis in the UTM projection;
- `shift_y`: when `shift` is `TRUE`, this will be the meters of the shift along y-axis in the UTM projection;
- `n_threads`: the maximum number of cores to be used in the smoothing and in the ISOFIT procedures;
- `aod_fixed`: when using the `isofit` step this will make you use the ISOFIT automatic estimation of Aerosol Optical Depth (AOD) or will try to fix the AOD in the atmospheric correction process to the best NASA Giovanni AOD value for each PRISMA image;
- `procedure_order`: a vector containing the list of all the processing steps in the order you want them. There is already a suggestion for L0, L1 and L2 products, but you can use it as LEGO building blocks.

After putting all these data inside then you can run everything clicking on `Source` in the upper right corner of the script.  
# DEFAULT USE CASES
Use cases for:
1) L0 products:
```r
   procedure_order <- c("inject","read","cloud","coreg","isofit","regrid","smooth","crop")
```
2) L1 products:
```r
   procedure_order <- c("read","cloud","coreg","isofit","regrid","smooth","crop")
```
3) L2 products:
```r
   procedure_order <- c("read","coreg","regrid","smooth","crop")
```
Just a note: PRISTAR-GEOSPEC will save each intermediate step, so don't worry if you want to stop before smoothing or cropping. 
# USE CASES
Use case:
1) you want read a L0 PRISMA image. Put `data_SubAcq3_C_SWIR_SURFACE-OBS_Part0_S11.h5` and `data_SubAcq3_C_VNIR_SURFACE-OBS_Part0_S11.he5` into `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. Then choose a L1 PRISMA image with same or similar view angle to the L0 product and put its `.he5` file inside the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. Then use PRISTAR-GEOSPEC with:
```r
   procedure_order <- c("inject","read")
```
2) you want to read a L1 or L2 PRISMA image in `.he5` format. Put your PRISMA file in `.he5` format in the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. Then use PRISTAR-GEOSPEC with:
```r
   procedure_order <- c("read")
```
3) you want to read, generate cloud mask and atmospherically correct an L1 PRISMA image. Put your PRISMA file in `.he5` format in the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. Put the dtm, the aspect and the slope rasters in the `config_folder/DEM` folder as separate files respectively with `dtm`, `aspect` and `slope` in their names. Then use PRISTAR-GEOSPEC with:
```r
   procedure_order <- c("read","cloud","isofit")
```
4) you want to read and coregister your L2 PRISMA image to a one-band Sentinel-2 image using gdalwarp. Put your PRISMA file in `.he5` format and Sentinel-2 one-band image in `.tif` format in the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. To optimize coregistration, choose the PRISMA band most suitable for the one-band you chose for the Sentinel-2 image (e.g. 52nd PRISMA image band if I chose the B8 Sentinel-2 band). Then use PRISTAR-GEOSPEC with:
```r
   procedure_order <- c("read","coreg")
   PRS_band_for_coreg <- 52
```
5) you want to read, generate cloud mask, coregister and atmospherically correct an L1 PRISMA image. Put your PRISMA file in `.he5` format and Sentinel-2 one-band image in `.tif` format in the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder. Put the dtm, the aspect and the slope rasters in the `config_folder/DEM` folder as separate files respectively with `dtm`, `aspect` and `slope` in their names. Choose the PRISMA band for coregistration and if you want a validation for coregistration. Then use PRISTAR-GEOSPEC with:
```r
   procedure_order <- c("read","cloud","coreg","isofit")
   PRS_band_for_coreg <- 52
```
# READING RESULTS
Inside the `config_folder/put_PRISMA_he5_and_S2_tif_here` folder you should have the list of PRISMA images folders (e.g L1 products):

<img width="680" height="141" alt="image" src="https://github.com/user-attachments/assets/f408c35f-edc4-4284-bec3-b26e0a878f71" />

If you open one of them you should find:

<img width="814" height="168" alt="image" src="https://github.com/user-attachments/assets/5ecf8b00-4a47-4752-b2d8-bbc12f8420a8" />

and if you open PRISTAR-processing folder you should find for the default L1 workflow:

<img width="497" height="167" alt="image" src="https://github.com/user-attachments/assets/93f3ae91-a749-44d5-bad1-caf40bba7f5c" />

Each folder contains the product in its name with the following naming convention:
- S: smoothing
- C: coregistration
- R: regridding
- O: orthoprojection
- T: trimming (cropping)
- A: atmospheric correction
- I: inject

that will appear in the filename in the order they have been performed from left to right.

# WHAT IF A PROBLEM
Open an Issue over the GitHub page & write email to the authors.
