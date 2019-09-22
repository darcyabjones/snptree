ARG IMAGE

FROM "${IMAGE}" as snpeff_builder

ARG SNPEFF_VERSION="4_3t"
ARG SNPEFF_URL="https://dl.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip"
ARG SNPEFF_PREFIX_ARG="/opt/snpeff/${SNPEFF_VERSION}"
ENV SNPEFF_PREFIX="${SNPEFF_PREFIX_ARG}"



WORKDIR /tmp
RUN  set -eux \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       unzip \
       wget \
       zlib1g-dev \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O snpeff.zip "${SNPEFF_URL}" \
  && unzip snpeff.zip \
  && mkdir -p "$(dirname ${SNPEFF_PREFIX})" \
  && mv snpEff "${SNPEFF_PREFIX}" \
  && sed -i 's~data.dir = ./data/~data.dir = /data/~' "${SNPEFF_PREFIX}/snpEff.config" \
  && rm -rf -- "${SNPEFF_PREFIX}/examples" "${SNPEFF_PREFIX}/galaxy" \
  && add_runtime_dep default-jre perl python

ARG SNPSIFT_SCRIPT
COPY "${SNPSIFT_SCRIPT}" "${SNPEFF_PREFIX}/scripts/SnpSift"

FROM "${IMAGE}"

ARG SNPEFF_VERSION="4_3t"
ARG SNPEFF_PREFIX_ARG="/opt/snpeff/${SNPEFF_VERSION}"
ENV SNPEFF_PREFIX="${SNPEFF_PREFIX_ARG}"
LABEL snpeff.version="${SNPEFF_VERSION}"

ENV PATH "${SNPEFF_PREFIX}/scripts:${PATH}"

COPY --from=snpeff_builder "${SNPEFF_PREFIX}" "${SNPEFF_PREFIX}"
COPY --from=snpeff_builder "${APT_REQUIREMENTS_FILE}" /build/apt/snpeff.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
