#
# Dockerfile
#
#

FROM rocker/rstudio:4.6.1
# https://rocker-project.org/images/versioned/rstudio.html

LABEL project="PRISTAR" \
      author="Emanuele Spirito" \
      maintainer="emanuelespirito@cnr.it" \
      image_name="PRISTAR" \
      version="0.4.0" \
      released="2026-07-01" \
      description="Rstudio"

USER root

# set graphical mode
ENV TERM=xterm
ENV DISPLAY=:1.0

# Set default locale
ENV LANG=C.UTF-8
ENV LANGUAGE=en_US:en

# Set default timezone
ENV TZ=UTC

# set environenmental variables
ENV EDITOR='nano'
ENV GDAL_NETCDF_CENTERLONG_180='NO'
ENV HDF5_USE_FILE_LOCKING=FALSE

# ARG DEBIAN_FRONTEND=noninteractive

# install all dependencies
RUN apt update -qq && DEBIAN_FRONTEND=noninteractive apt install -yq \ 
    build-essential \ 
    gfortran \
    libabsl-dev \ 
    cmake \ 
    libgdal-dev \ 
    gdal-bin \ 
    libgeos-dev \ 
    libproj-dev \ 
    libsqlite3-dev \ 
    libudunits2-dev \ 
    libicu-dev \ 
    libuv1-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

#Base packages
RUN R --vanilla -e 'install.packages(c("pak"), repos = "https://cloud.r-project.org")'

#Main packages and Github packages
RUN R --vanilla -e 'pak::pkg_install(c("here", "remotes", "httr2", "optparse", "tidytable", "data.table", "lubridate", "hms", \
                                       "sf", "terra", "exactextractr", "mmand", "tidyterra", "satellite", \
                                       "irea-cnr-mi/prismaread", "kamapu/spatialist")); pak::pak_cleanup(force = TRUE)'

# Disable automatic Ubuntu updates
RUN echo 'APT::Periodic::Update-Package-Lists "0";' > /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::Download-Upgradeable-Packages "0";' >> /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::AutocleanInterval "0";' >> /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::Unattended-Upgrade "0";' >> /etc/apt/apt.conf.d/20auto-upgrades



CMD ["/init"]
