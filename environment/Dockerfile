FROM registry.codeocean.com/codeocean/julia:1.12.1-22.04.1

ARG DEBIAN_FRONTEND=noninteractive

COPY Project.toml /julia-env/Project.toml
COPY postInstall /
RUN /postInstall
