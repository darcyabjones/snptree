ARG IMAGE

FROM "${IMAGE}" as nwutils_builder

ARG NWUTILS_VERSION
ARG NWUTILS_URL
ARG NWUTILS_PREFIX_ARG
ENV NWUTILS_PREFIX="${NWUTILS_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y \
       bison \
       build-essential \
       ca-certificates \
       file \
       flex-old \
       guile-2.0-dev \
       libxml2-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O nwutils.tar.gz "${NWUTILS_URL}" \
  && tar -zxf nwutils.tar.gz \
  && cd newick-utils-*/ \ 
  && ./configure --prefix="${NWUTILS_PREFIX}" \
  && make \
  && make check \
  && make install \
  && add_runtime_dep guile-2.0 libxml2


FROM "${IMAGE}"

ARG NWUTILS_VERSION
ARG NWUTILS_PREFIX_ARG
ENV NWUTILS_PREFIX="${NWUTILS_PREFIX_ARG}"

LABEL nwutils.version="${NWUTILS_VERSION}"

ENV PATH="${NWUTILS_PREFIX}/bin:${PATH}"

COPY --from=nwutils_builder "${NWUTILS_PREFIX}" "${NWUTILS_PREFIX}"
COPY --from=nwutils_builder "${APT_REQUIREMENTS_FILE}" /build/apt/nwutils.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
