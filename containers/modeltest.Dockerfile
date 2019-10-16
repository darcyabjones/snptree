ARG IMAGE

FROM "${IMAGE}" as modeltest_builder

ARG MODELTEST_COMMIT
ARG MODELTEST_REPO
ARG MODELTEST_PREFIX_ARG="/opt/modeltest/${MODELTEST_COMMIT}"
ENV MODELTEST_PREFIX="${MODELTEST_PREFIX_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       autogen \
       autoconf \
       automake \
       autotools-dev \
       bison \
       build-essential \
       ca-certificates \
       cmake \
       flex \
       git \
       libgtest-dev \
       libmpich-dev \
       libtool \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone --recursive "${MODELTEST_REPO}" . \
  && mkdir build \
  && cd build \
  && cmake \
       -DUSE_MPI=ON \
       -DUSE_PTHREADS=ON \
       -DUSE_GUI=OFF \
       -DSTATIC_BUILD=ON \
       -DENABLE_MODELTEST_SIMD=ON \
       -DENABLE_PLLMOD_SIMD=ON \
       -DCMAKE_INSTALL_PREFIX="${MODELTEST_PREFIX}" \
       .. \
  && make \
  && mkdir -p "${MODELTEST_PREFIX}/bin" \
  && cp ../bin/modeltest-ng-mpi "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi-avx" \
  && make clean \
  && cmake \
       -DUSE_MPI=ON \
       -DUSE_PTHREADS=ON \
       -DUSE_GUI=OFF \
       -DSTATIC_BUILD=ON \
       -DENABLE_MODELTEST_SIMD=OFF \
       -DENABLE_PLLMOD_SIMD=OFF \
       -DENABLE_SSE=false \
       -DENABLE_AVX=false \
       -DENABLE_AVX2=false \
       -DCMAKE_INSTALL_PREFIX="${MODELTEST_PREFIX}" \
       .. \
  && make \
  && mkdir -p "${MODELTEST_PREFIX}/bin" \
  && cp ../bin/modeltest-ng-mpi "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi-noavx" \
  && make clean \
  && cmake \
       -DUSE_MPI=OFF \
       -DUSE_PTHREADS=ON \
       -DUSE_GUI=OFF \
       -DSTATIC_BUILD=ON \
       -DENABLE_MODELTEST_SIMD=ON \
       -DENABLE_PLLMOD_SIMD=ON \
       -DCMAKE_INSTALL_PREFIX="${MODELTEST_PREFIX}" \
       .. \
  && make \
  && cp ../bin/modeltest-ng-static "${MODELTEST_PREFIX}/bin/modeltest-ng-avx" \
  && make clean \
  && cmake \
       -DUSE_MPI=OFF \
       -DUSE_PTHREADS=ON \
       -DUSE_GUI=OFF \
       -DSTATIC_BUILD=ON \
       -DENABLE_MODELTEST_SIMD=OFF \
       -DENABLE_PLLMOD_SIMD=OFF \
       -DCMAKE_INSTALL_PREFIX="${MODELTEST_PREFIX}" \
       .. \
  && make \
  && mkdir -p "${MODELTEST_PREFIX}/bin" \
  && cp ../bin/modeltest-ng-static "${MODELTEST_PREFIX}/bin/modeltest-ng-noavx" \
  && make clean \
  && echo '#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${MODELTEST_PREFIX}/bin/modeltest-ng-avx" "$@" \n\
else \n\
    exec "${MODELTEST_PREFIX}/bin/modeltest-ng-noavx" "$@" \n\
fi \n\
' > "${MODELTEST_PREFIX}/bin/modeltest-ng" \
  && chmod a+x "${MODELTEST_PREFIX}/bin/modeltest-ng" \
  && echo '#!/usr/bin/env bash\n\
if $(grep -q -E '^flags.+avx2' /proc/cpuinfo); then \n\
    exec "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi-avx" "$@" \n\
else \n\
    exec "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi-noavx" "$@" \n\
fi \n\
' > "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi" \
  && chmod a+x "${MODELTEST_PREFIX}/bin/modeltest-ng-mpi" \
  && add_runtime_dep mpich


FROM "${IMAGE}"

ARG MODELTEST_COMMIT
ARG MODELTEST_PREFIX_ARG="/opt/modeltest/${MODELTEST_COMMIT}"
ENV MODELTEST_PREFIX="${MODELTEST_PREFIX_ARG}"

LABEL modeltest.version="${MODELTEST_COMMIT}"

ENV PATH "${MODELTEST_PREFIX}/bin:${PATH}"

COPY --from=modeltest_builder "${MODELTEST_PREFIX}" "${MODELTEST_PREFIX}"
COPY --from=modeltest_builder "${APT_REQUIREMENTS_FILE}" /build/apt/modeltest.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
