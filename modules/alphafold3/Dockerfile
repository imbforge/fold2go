FROM ghcr.io/prefix-dev/pixi:jammy

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV DEBIAN_FRONTEND="noninteractive" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" XLA_FLAGS="--xla_gpu_enable_triton_gemm=false" XLA_PYTHON_CLIENT_MEM_FRACTION=0.95 XLA_PYTHON_CLIENT_PREALLOCATE=true

RUN apt-get update && apt-get install -y --no-install-recommends build-essential ca-certificates git zlib1g-dev

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN bash -c "install -D <(pixi shell-hook) /.singularity.d/env/99-pixi.sh" && pixi clean cache --yes

RUN /app/alphafold/.pixi/envs/default/bin/build_data
