FROM ghcr.io/prefix-dev/pixi:jammy AS build

ENV DEBIAN_FRONTEND="noninteractive" CONDA_OVERRIDE_CUDA="12.0"

WORKDIR /app/openfold

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

FROM ghcr.io/prefix-dev/pixi:jammy AS openfold

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute" KMP_AFFINITY="none" TORCH_CUDA_ARCH_LIST="8.0;8.6;9.0"

COPY --from=build /app/openfold /app/openfold
COPY --from=build /app/openfold/99-pixi.sh /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["run_openfold", "predict"]
