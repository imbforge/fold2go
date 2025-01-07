FROM ghcr.io/prefix-dev/pixi:jammy AS build

ENV DEBIAN_FRONTEND="noninteractive"

RUN apt-get update && apt-get install -y --no-install-recommends build-essential ca-certificates git

WORKDIR /app/boltz

COPY pixi.toml pixi.lock .

RUN printf '#!/bin/sh\n%s\nexec "$@"' "$(pixi shell-hook)" > 99-pixi.sh && chmod +x 99-pixi.sh

FROM docker.io/nvidia/cuda:12.6.0-base-ubuntu22.04 AS boltz

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

ENV NVIDIA_VISIBLE_DEVICES="all" NVIDIA_DRIVER_CAPABILITIES="compute"

COPY --from=build /app/boltz /app/boltz
COPY --from=build /app/boltz/99-pixi.sh /.singularity.d/env/99-pixi.sh

ENTRYPOINT ["/.singularity.d/env/99-pixi.sh"]

CMD ["boltz", "predict"]