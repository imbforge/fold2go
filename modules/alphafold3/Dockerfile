FROM ghcr.io/prefix-dev/pixi:jammy AS build

ENV DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends build-essential ca-certificates git zlib1g-dev

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

RUN /app/alphafold/.pixi/envs/default/bin/build_data

FROM docker.io/nvidia/cuda:12.6.0-base-ubuntu22.04 AS alphafold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV PATH="/app/alphafold:$PATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" XLA_FLAGS="--xla_gpu_enable_triton_gemm=false" XLA_PYTHON_CLIENT_PREALLOCATE=false TF_FORCE_UNIFIED_MEMORY=true XLA_PYTHON_CLIENT_MEM_FRACTION=3.2

COPY --from=build /app/alphafold /app/alphafold

WORKDIR /app/alphafold

COPY --from=build /app/alphafold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ADD --chmod=0755 https://raw.githubusercontent.com/google-deepmind/alphafold3/e7510b111740b7e97b25310c531081fdcdd3cdf7/run_alphafold.py .

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["python", "/app/alphafold/run_alphafold.py"]