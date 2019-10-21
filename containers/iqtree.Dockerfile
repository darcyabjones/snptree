ARG IMAGE

FROM "${IMAGE}" as iqtree_builder

ARG IQTREE_VERSION
ARG IQTREE_URL
ARG IQTREE_PREFIX_ARG="/opt/iqtree/${IQTREE_VERSION}"
ENV IQTREE_PREFIX="${IQTREE_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       cmake \
       libboost-dev \
       libeigen3-dev \
       libmpich-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && wget -O iqtree.tar.gz "${IQTREE_URL}" \
  && tar zxf iqtree.tar.gz \
  && cd IQ-TREE*/ \
  && mkdir build \
  && cd build \
  && cmake -DIQTREE_FLAGS="omp novx" -DCMAKE_INSTALL_PREFIX="${IQTREE_PREFIX}" .. \
  && make \
  && make install \
  && mv "${IQTREE_PREFIX}/bin/iqtree" "${IQTREE_PREFIX}/bin/iqtree-sse" \
  && cmake -DIQTREE_FLAGS="omp-mpi novx" -DCMAKE_INSTALL_PREFIX="${IQTREE_PREFIX}" .. \
  && make \
  && make install \
  && mv "${IQTREE_PREFIX}/bin/iqtree-mpi" "${IQTREE_PREFIX}/bin/iqtree-sse-mpi" \
  && cmake -DIQTREE_FLAGS="omp-mpi" -DCMAKE_INSTALL_PREFIX="${IQTREE_PREFIX}" .. \
  && make \
  && make install \
  && mv "${IQTREE_PREFIX}/bin/iqtree-mpi" "${IQTREE_PREFIX}/bin/iqtree-avx-mpi" \
  && cmake -DIQTREE_FLAGS="omp" -DCMAKE_INSTALL_PREFIX="${IQTREE_PREFIX}" .. \
  && make \
  && make install \
  && mv "${IQTREE_PREFIX}/bin/iqtree" "${IQTREE_PREFIX}/bin/iqtree-avx" \
  && echo '#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${IQTREE_PREFIX}/bin/iqtree-avx" "$@" \n\
else \n\
    exec "${IQTREE_PREFIX}/bin/iqtree-sse" "$@" \n\
fi \n\
' > "${IQTREE_PREFIX}/bin/iqtree" \
  && chmod a+x "${IQTREE_PREFIX}/bin/iqtree" \
  && echo $'#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${IQTREE_PREFIX}/bin/iqtree-avx-mpi" "$@" \n\
else \
    exec "${IQTREE_PREFIX}/bin/iqtree-sse-mpi" "$@" \n\
fi \n\
' > "${IQTREE_PREFIX}/bin/iqtree-mpi" \
  && chmod a+x "${IQTREE_PREFIX}/bin/iqtree-mpi" \
  && add_runtime_dep libgomp1 mpich

# CA cert stuff sometime required for git clone https


FROM "${IMAGE}"

ARG IQTREE_VERSION
ARG IQTREE_PREFIX_ARG="/opt/iqtree/${IQTREE_VERSION}"
ENV IQTREE_PREFIX="${IQTREE_PREFIX_ARG}"
LABEL iqtree.version="${IQTREE_VERSION}"

ENV PATH "${IQTREE_PREFIX}/bin:${PATH}"

COPY --from=iqtree_builder "${IQTREE_PREFIX}" "${IQTREE_PREFIX}"
COPY --from=iqtree_builder "${APT_REQUIREMENTS_FILE}" /build/apt/iqtree.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
