FROM ghcr.io/prefix-dev/pixi:jammy as build

ENV DEBIAN_FRONTEND="noninteractive" CONDA_OVERRIDE_CUDA="12.2.2"

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

ADD https://github.com/google-deepmind/alphafold/archive/c77e5d2a8961d1a353632c462914ff0a32a950f6.tar.gz /app/alphafold/alphafold.tar.gz
ADD https://github.com/sokrypton/ColabFold/raw/v1.6.3/colabfold/openstructure/stereo_chemical_props.txt /app/alphafold/alphafold/common/

ADD http://eddylab.org/software/hmmer/hmmer-3.4.tar.gz /tmp/

# add assorted easel and jackhmmer patches from alphafold3 and easel repositories
ADD https://github.com/EddyRivasLab/easel/commit/c0497d24c043c90129aa8f812f4bc8d7aa89e5c5.patch /tmp/
ADD https://github.com/google-deepmind/alphafold3/raw/v3.0.5/docker/easel_fast_fasta_parsing.patch /tmp/
ADD https://github.com/google-deepmind/alphafold3/raw/v3.0.5/docker/jackhmmer_parallel_fasta.patch /tmp/
ADD https://github.com/google-deepmind/alphafold3/raw/v3.0.5/docker/jackhmmer_seq_limit.patch /tmp/

RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
    && tar xzf alphafold.tar.gz --strip-components=1 \
    && chmod +x run_alphafold.py \
    && rm alphafold.tar.gz

WORKDIR /tmp

RUN tar xzf hmmer-3.4.tar.gz \
    && patch hmmer-3.4/src/jackhmmer.c /tmp/jackhmmer_seq_limit.patch \
    && patch hmmer-3.4/src/jackhmmer.c /tmp/jackhmmer_parallel_fasta.patch \
    && patch hmmer-3.4/easel/esl_sqio_ascii.c /tmp/easel_fast_fasta_parsing.patch \
    && patch hmmer-3.4/easel/esl_workqueue.c /tmp/c0497d24c043c90129aa8f812f4bc8d7aa89e5c5.patch

WORKDIR /tmp/hmmer-3.4

RUN ./configure --prefix=/app/alphafold/.pixi/envs/default && make install && cd easel && make install

FROM docker.io/nvidia/cuda:12.2.2-base-ubuntu22.04 AS alphafold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV PATH="/opt/conda/bin:/app/alphafold:$PATH" PYTHONPATH="/app/alphafold:$PYTHONPATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

COPY --from=build /app/alphafold /app/alphafold
COPY --from=build /app/alphafold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["python", "/app/alphafold/run_alphafold.py"]
