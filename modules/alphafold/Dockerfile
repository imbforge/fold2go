FROM ghcr.io/prefix-dev/pixi:jammy-cuda-12.2.2

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV DEBIAN_FRONTEND="noninteractive" PATH="/opt/conda/bin:/app/alphafold:$PATH" PYTHONPATH="/app/alphafold:$PYTHONPATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN bash -c "install -D <(pixi shell-hook) /.singularity.d/env/99-pixi.sh" && pixi clean cache --yes

ADD https://github.com/google-deepmind/alphafold/archive/f251de6613cb478207c732bf9627b1e853c99c2f.tar.gz /app/alphafold/alphafold.tar.gz
RUN tar xvz --strip-components=1 alphafold.tar.gz && chmod +x run_alphafold.py && rm alphafold.tar.gz
ADD --chmod=664 https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt /app/alphafold/alphafold/common/
