FROM ghcr.io/prefix-dev/pixi:jammy

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV DEBIAN_FRONTEND="noninteractive" PATH="/opt/conda/bin:$PATH" NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=4.0 OPENMM_CPU_THREADS=8

COPY pixi.toml pixi.lock .

RUN mkdir -p /.singularity.d/env && printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > /.singularity.d/env/99-pixi.sh && chmod +x /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["colabfold_batch"]