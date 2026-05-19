FROM ghcr.io/prefix-dev/pixi:jammy as build

ENV DEBIAN_FRONTEND="noninteractive" CONDA_OVERRIDE_CUDA="12.2.2"

WORKDIR /app/alphafold

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

ADD https://github.com/google-deepmind/alphafold/archive/f251de6613cb478207c732bf9627b1e853c99c2f.tar.gz /app/alphafold/alphafold.tar.gz
RUN tar -xvzf alphafold.tar.gz --strip-components=1 && chmod +x run_alphafold.py && rm alphafold.tar.gz
ADD --chmod=664 https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt /app/alphafold/alphafold/common/

FROM docker.io/nvidia/cuda:12.2.2-base-ubuntu22.04 AS alphafold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV PATH="/opt/conda/bin:/app/alphafold:$PATH" PYTHONPATH="/app/alphafold:$PYTHONPATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

COPY --from=build /app/alphafold /app/alphafold
COPY --from=build /app/alphafold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["python", "/app/alphafold/run_alphafold.py"]
