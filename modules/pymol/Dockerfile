FROM ghcr.io/prefix-dev/pixi:jammy

LABEL org.opencontainers.image.authors="huether@imb-mainz.de"

WORKDIR /opt

COPY pixi.toml pixi.lock .

RUN bash -c "install -D <(pixi shell-hook) /.singularity.d/env/99-pixi.sh" && pixi clean cache --yes
