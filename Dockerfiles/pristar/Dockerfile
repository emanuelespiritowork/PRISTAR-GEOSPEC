#
# Dockerfile
#
#

FROM rocker/rstudio:4.5.1
# https://rocker-project.org/images/versioned/rstudio.html

LABEL project="PRISTAR" \
      author="Emanuele Spirito" \
      maintainer="XXXX@gmail.com" \
      image_name="PRISTAR" \
      version="test" \
      released="2026-06-30" \
      software_versions="Ubuntu 24.04 XXXXXXXXXXXXX" \
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
RUN R --vanilla -e 'install.packages(c("remotes", "pak"), repos = "https://cloud.r-project.org", Ncpus = 8)'

#Main packages
RUN R --vanilla -e 'pak::pkg_install(c("here", "httr2", "optparse@1.7.5", "tidytable@0.11.2", "data.table@1.17.6", "sf@1.0-24", "terra@1.8-54", "exactextractr@0.10.0", "mmand@1.6.3","tidyterra@0.7.2")); pak::pak_cleanup(force = TRUE)'

#Github packages
RUN R --vanilla -e 'pak::pkg_install(c("irea-cnr-mi/prismaread", "kamapu/spatialist")); pak::pak_cleanup(force = TRUE)'

# Disable automatic Ubuntu updates
RUN echo 'APT::Periodic::Update-Package-Lists "0";' > /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::Download-Upgradeable-Packages "0";' >> /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::AutocleanInterval "0";' >> /etc/apt/apt.conf.d/20auto-upgrades && echo 'APT::Periodic::Unattended-Upgrade "0";' >> /etc/apt/apt.conf.d/20auto-upgrades



CMD ["/init"]
