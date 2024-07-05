FROM ghcr.io/prefix-dev/pixi:jammy

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

RUN apt-get update && apt-get install --no-install-recommends -y libgl1

WORKDIR /opt

COPY pixi.toml pixi.lock .

RUN bash -c "install -D <(pixi shell-hook) /.singularity.d/env/99-pixi.sh" && pixi clean cache --yes
