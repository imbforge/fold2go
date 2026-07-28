ARG OPENFOLD3_VERSION=0.4.4

FROM ghcr.io/prefix-dev/pixi:jammy AS build

ARG OPENFOLD3_VERSION

ENV DEBIAN_FRONTEND="noninteractive" CONDA_OVERRIDE_CUDA="12.0"

WORKDIR /app/openfold

ADD https://github.com/aqlaboratory/openfold-3/archive/refs/tags/${OPENFOLD3_VERSION}.tar.gz /app/openfold/

RUN tar xzf ${OPENFOLD3_VERSION}.tar.gz --strip-components=1 && printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook -e openfold3-cuda12)" > 99-pixi.sh && chmod +x 99-pixi.sh

FROM ghcr.io/prefix-dev/pixi:jammy AS openfold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" KMP_AFFINITY="none" TORCH_CUDA_ARCH_LIST="8.0;8.6;9.0"

COPY --from=build /app/openfold/.pixi /app/openfold/pixi.toml /app/openfold/pixi.lock /app/openfold/
COPY --from=build /app/openfold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["run_openfold", "predict"]
