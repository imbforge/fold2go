FROM ghcr.io/prefix-dev/pixi:jammy AS build

ENV DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends build-essential ca-certificates git zlib1g-dev

# create pixi environment for alphafold3

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

RUN pixi run build_data

# build hmmer with af3 seq_limit patch to reduce peak RAM usage

WORKDIR /tmp

ADD https://raw.githubusercontent.com/google-deepmind/alphafold3/751a4b8612d0d53de8f6e1830c8f726e873a55cf/docker/jackhmmer_seq_limit.patch /tmp/
ADD http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz /tmp/

RUN tar xzf hmmer-3.4.tar.gz && patch -p0 < jackhmmer_seq_limit.patch

WORKDIR /tmp/hmmer-3.4

RUN ./configure --prefix=/app/alphafold/.pixi/envs/default && make install && cd easel && make install

FROM docker.io/nvidia/cuda:12.6.3-base-ubuntu22.04 AS alphafold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV PATH="/app/alphafold:$PATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" XLA_FLAGS="--xla_gpu_enable_triton_gemm=false" XLA_PYTHON_CLIENT_PREALLOCATE=false TF_FORCE_UNIFIED_MEMORY=true XLA_PYTHON_CLIENT_MEM_FRACTION=3.2

COPY --from=build /app/alphafold /app/alphafold

WORKDIR /app/alphafold

COPY --from=build /app/alphafold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ADD --chmod=0755 https://raw.githubusercontent.com/google-deepmind/alphafold3/f6032f0f6ad78a93ad405675d2c8fa27f3cec093/run_alphafold.py .

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["python", "/app/alphafold/run_alphafold.py"]
