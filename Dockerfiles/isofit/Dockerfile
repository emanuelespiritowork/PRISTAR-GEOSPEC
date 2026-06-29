FROM --platform=$BUILDPLATFORM python:3.12-slim
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

USER root
RUN apt-get update && \
    apt-get install -y --no-install-suggests --no-install-recommends \
      g++ \
      gfortran \
      make \
      libgsl-dev \
      libnetcdf-dev && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root

COPY . ISOFIT/
RUN cd ISOFIT && \
    uv sync && \
    uv run isofit -b /root download all && \
    uv run isofit build && \
    echo "source /root/ISOFIT/.venv/bin/activate" >> ~/.bashrc

# Add FastAPI to the existing venv
RUN cd ISOFIT && uv pip install fastapi "uvicorn[standard]"

COPY api.py /root/api.py

ENV SHELL="/bin/bash"
ENV PATH="/root/ISOFIT/.venv/bin:$PATH"
ENV VIRTUAL_ENV="/root/ISOFIT/.venv"

# Ray Dashboard port
EXPOSE 8265
# API port
EXPOSE 8000

CMD ["uvicorn", "api:app", "--host", "0.0.0.0", "--port", "8000", "--app-dir", "/root"]