ARG IMAGE

FROM "${IMAGE}" as raxml_builder

ARG RAXML_COMMIT
ARG RAXML_REPO
ARG RAXML_PREFIX_ARG="/opt/raxml/${RAXML_COMMIT}"
ENV RAXML_PREFIX="${RAXML_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       bison \
       build-essential \
       ca-certificates \
       cmake \
       flex \
       git \
       libgtest-dev \
       libgmp3-dev \
       libmpich-dev \
       wget \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone --recursive "${RAXML_REPO}" . \
  && mkdir build \
  && cd build \
  && cmake \
       -DUSE_MPI=ON \
       -DUSE_PTHREADS=ON \
       -DUSE_GMP=ON \
       -DSTATIC_BUILD=ON \
       -DENABLE_RAXML_SIMD=ON \
       -DENABLE_PLLMOD_SIMD=ON \
       -DCMAKE_INSTALL_PREFIX="${RAXML_PREFIX}" \
       .. \
  && make \
  && make install \
  && mv "${RAXML_PREFIX}/bin/raxml-ng-mpi" "${RAXML_PREFIX}/bin/raxml-ng-mpi-avx" \
  && make clean \
  && cmake \
       -DUSE_MPI=ON \
       -DUSE_PTHREADS=ON \
       -DUSE_GMP=ON \
       -DSTATIC_BUILD=ON \
       -DENABLE_RAXML_SIMD=OFF \
       -DENABLE_PLLMOD_SIMD=OFF \
       -DENABLE_SSE=false \
       -DENABLE_AVX=false \
       -DENABLE_AVX2=false \
       -DCMAKE_INSTALL_PREFIX="${RAXML_PREFIX}" \
       .. \
  && make \
  && make install \
  && mv "${RAXML_PREFIX}/bin/raxml-ng-mpi" "${RAXML_PREFIX}/bin/raxml-ng-mpi-noavx" \
  && cmake \
       -DUSE_MPI=OFF \
       -DUSE_PTHREADS=ON \
       -DUSE_GMP=ON \
       -DSTATIC_BUILD=ON \
       -DENABLE_RAXML_SIMD=OFF \
       -DENABLE_PLLMOD_SIMD=OFF \
       -DENABLE_SSE=false \
       -DENABLE_AVX=false \
       -DENABLE_AVX2=false \
       -DCMAKE_INSTALL_PREFIX="${RAXML_PREFIX}" \
       .. \
  && make \
  && make install \
  && mv "${RAXML_PREFIX}/bin/raxml-ng-static" "${RAXML_PREFIX}/bin/raxml-ng-noavx" \
  && echo '#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${RAXML_PREFIX}/bin/raxml-ng-avx" "$@" \n\
else \n\
    exec "${RAXML_PREFIX}/bin/raxml-ng-noavx" "$@" \n\
fi \n\
' > "${RAXML_PREFIX}/bin/raxml-ng" \
  && chmod a+x "${RAXML_PREFIX}/bin/raxml-ng" \
  && echo '#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${RAXML_PREFIX}/bin/raxml-ng-mpi-avx" "$@" \n\
else \n\
    exec "${RAXML_PREFIX}/bin/raxml-ng-mpi-noavx" "$@" \n\
fi \n\
' > "${RAXML_PREFIX}/bin/raxml-ng-mpi" \
  && chmod a+x "${RAXML_PREFIX}/bin/raxml-ng-mpi" \
  && add_runtime_dep mpich libgmp10 libgmpxx4ldbl


FROM "${IMAGE}"

ARG RAXML_COMMIT
ARG RAXML_PREFIX_ARG="/opt/raxml/${RAXML_COMMIT}"
ENV RAXML_PREFIX="${RAXML_PREFIX_ARG}"

LABEL raxml.version="${RAXML_COMMIT}"

ENV PATH "${RAXML_PREFIX}/bin:${PATH}"

COPY --from=raxml_builder "${RAXML_PREFIX}" "${RAXML_PREFIX}"
COPY --from=raxml_builder "${APT_REQUIREMENTS_FILE}" /build/apt/raxml.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
