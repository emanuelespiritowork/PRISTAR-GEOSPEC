#
# PRISTAR — R + ISOFIT + AROSICS in a single image
#
# Base: rocker/rstudio (Ubuntu noble). For the CLI-only version later,
# switch to rocker/r-ver:4.6.1 and change CMD (see bottom).
#
FROM rocker/rstudio:4.6.1

LABEL project="PRISTAR" \
      author="Emanuele Spirito" \
      maintainer="emanuelespirito@cnr.it" \
      image_name="PRISTAR" \
      version="1.0.0" \
      description="R + ISOFIT + AROSICS (isolated uv venvs)"

USER root

ENV TERM=xterm \
    DISPLAY=:1.0 \
    LANG=C.UTF-8 \
    LANGUAGE=en_US:en \
    TZ=UTC \
    EDITOR=nano \
    GDAL_NETCDF_CENTERLONG_180=NO \
    HDF5_USE_FILE_LOCKING=FALSE

# ---------------------------------------------------------------------------
# 1. System dependencies (union of the three images)
#    R/geo:   gdal, geos, proj, udunits, ...
#    ISOFIT:  gfortran, gsl, netcdf (+ build-essential for g++/make)
#    AROSICS: python3-dev (GDAL bindings build), libspatialindex
# ---------------------------------------------------------------------------
RUN apt-get update -qq && DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        build-essential \
        gfortran \
        cmake \
        libabsl-dev \
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
        libgsl-dev \
        libnetcdf-dev \
        libspatialindex-dev \
        python3 \
        python3-dev \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. R packages
# ---------------------------------------------------------------------------
RUN R -e 'install.packages("pak")'

RUN R -e 'pak::pkg_install(c( \
        "here", "remotes", "httr2", "optparse", "processx", \
        "tidytable", "data.table", "lubridate", "hms", \
        "sf", "terra", "exactextractr", "mmand", "tidyterra", "satellite", \
        "kamapu/spatialist", "emanuelespiritowork/prismaread")); \
      pak::pak_cleanup(force = TRUE)'

# ---------------------------------------------------------------------------
# 3. uv (pinned). Everything under /opt, NOT /root, so the runtime user
#    (rstudio / PUID) can read it. Venvs are NOT added to PATH.
# ---------------------------------------------------------------------------
COPY --from=ghcr.io/astral-sh/uv:0.11.24 /uv /uvx /bin/

ENV UV_PYTHON_INSTALL_DIR=/opt/uv/python \
    UV_LINK_MODE=copy \
    UV_NO_CACHE=1

# ---------------------------------------------------------------------------
# 4. ISOFIT venv  ->  /opt/isofit/.venv
#    Source stays in /opt/isofit/src (uv sync installs it as editable).
#    Data (sRTMnet, 6S, examples) is baked into the image at /opt/isofit/data
#    so the image is complete. Do NOT mount a volume on /opt/isofit/data.
#    Always call isofit with `-b /opt/isofit/data` at runtime.
# ---------------------------------------------------------------------------
COPY ISOFIT/ /opt/isofit/src/
RUN cd /opt/isofit/src \
    && UV_PROJECT_ENVIRONMENT=/opt/isofit/.venv uv sync --python /usr/bin/python3

ENV ISOFIT_DATA=/opt/isofit/data

# Download only what the pipeline needs:
#   essentials = data + sixs (downloaded AND compiled here) + srtmnet
#   surface    = EMIT surface prior files
# Skipped: examples, plots, libradtran. `isofit build` only builds the
# examples, so it is not needed anymore.
# `validate essentials` only prints a report in the build log (exit code 0).
RUN /opt/isofit/.venv/bin/isofit -b ${ISOFIT_DATA} download essentials \
    && /opt/isofit/.venv/bin/isofit -b ${ISOFIT_DATA} download surface \
    && /opt/isofit/.venv/bin/isofit -b ${ISOFIT_DATA} validate essentials \
    && chmod -R a+rX ${ISOFIT_DATA}

# ---------------------------------------------------------------------------
# 5. AROSICS venv  ->  /opt/arosics/.venv
#    GDAL Python bindings are compiled inside the venv against the system
#    libgdal (same one terra/sf use), instead of --system-site-packages.
#    numpy<2: GDAL < 3.9 bindings do not build against numpy 2.
# ---------------------------------------------------------------------------
RUN uv venv --python /usr/bin/python3 /opt/arosics/.venv \
    && echo "numpy<2" > /opt/arosics/constraints.txt \
    && uv pip install --python /opt/arosics/.venv/bin/python -c /opt/arosics/constraints.txt \
        numpy setuptools wheel \
    && uv pip install --python /opt/arosics/.venv/bin/python -c /opt/arosics/constraints.txt \
        --no-build-isolation "gdal==$(gdal-config --version)" \
    && uv pip install --python /opt/arosics/.venv/bin/python -c /opt/arosics/constraints.txt \
        cartopy \
        geopandas \
        "joblib>=1.3.0" \
        matplotlib \
        pandas \
        "pyproj>2.2.0" \
        "scikit-image>=0.21.0" \
        shapely \
        h5py \
        "arosics==1.13.2" \
    && /opt/arosics/.venv/bin/python -c "from osgeo import gdal, gdal_array; import arosics; print('GDAL', gdal.__version__, '| arosics', arosics.__version__)"

# ---------------------------------------------------------------------------
# 6. Launchers: `isofit-env <cmd> [args]` / `arosics-env <cmd> [args]`
#    Activate the venv for that single process only, then exec the command.
# ---------------------------------------------------------------------------
RUN for env in isofit arosics; do \
      printf '#!/bin/sh\nexport VIRTUAL_ENV=/opt/%s/.venv\nexport PATH="$VIRTUAL_ENV/bin:$PATH"\nexec "$@"\n' "$env" \
        > /usr/local/bin/${env}-env; \
      chmod 755 /usr/local/bin/${env}-env; \
    done

# Disable automatic Ubuntu updates
RUN printf '%s\n' \
      'APT::Periodic::Update-Package-Lists "0";' \
      'APT::Periodic::Download-Upgradeable-Packages "0";' \
      'APT::Periodic::AutocleanInterval "0";' \
      'APT::Periodic::Unattended-Upgrade "0";' \
      > /etc/apt/apt.conf.d/20auto-upgrades

# rocker/rstudio: /init starts RStudio Server.
# With rocker/r-ver as base: CMD ["R"] (or your Rscript entrypoint).
CMD ["/init"]