FROM ghcr.io/osgeo/gdal:ubuntu-small-3.13.1

# libspatialindex (per rtree/geopandas) e build tools
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3-dev \
        build-essential \
        libspatialindex-dev \
    && rm -rf /var/lib/apt/lists/*

COPY --from=ghcr.io/astral-sh/uv:0.11.24 /uv /uvx /bin/

ENV VIRTUAL_ENV="/opt/arosics/.venv"
RUN uv venv --system-site-packages $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN uv pip install \
        cartopy \
        geopandas \
        "joblib>=1.3.0" \
        matplotlib \
        numpy \
        pandas \
        "pyproj>2.2.0" \
        "scikit-image>=0.21.0" \
        shapely \
        h5py \
        "arosics>=1.3.0"

RUN uv pip install fastapi "uvicorn[standard]"

COPY api.py /app/api.py
WORKDIR /app
EXPOSE 8000

CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8000"]