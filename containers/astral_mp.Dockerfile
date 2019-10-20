ARG IMAGE

FROM "${IMAGE}" as astral_builder

ARG ASTRAL_COMMIT
ARG ASTRAL_REPO
ARG ASTRAL_PREFIX_ARG
ARG ASTRAL_JAR_ARG
ENV ASTRAL_PREFIX="${ASTRAL_PREFIX_ARG}"
ENV ASTRAL_JAR="${ASTRAL_JAR_ARG}"

WORKDIR /tmp
RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt-get install -y --no-install-recommends \
       build-essential \
       ca-certificates \
       default-jdk-headless \
       git \
       unzip \
       zip \
  && rm -rf /var/lib/apt/lists/* \
  && update-ca-certificates \
  && git clone "${ASTRAL_REPO}" astral \
  && cd astral \
  && git checkout "${ASTRAL_COMMIT}" \
  && rm Astral.*.zip \
  && ./make.sh \
  && unzip -o Astral.*.zip \
  && mkdir -p "${ASTRAL_PREFIX}" \
  && mv Astral/* "${ASTRAL_PREFIX}" \
  && chmod a+rx "${ASTRAL_JAR}" \
  && rm -rf -- "${ASTRAL_PREFIX}/test_data" "${ASTRAL_PREFIX}"/*.pdf \
  && add_runtime_dep default-jre-headless libjsap-java


FROM "${IMAGE}"

ARG ASTRAL_COMMIT
ARG ASTRAL_PREFIX_ARG
ARG ASTRAL_JAR_ARG
ENV ASTRAL_PREFIX="${ASTRAL_PREFIX_ARG}"
ENV ASTRAL_JAR="${ASTRAL_JAR_ARG}"

ENV LD_LIBRARY_PATH="${ASTRAL_PREFIX}/lib:${LD_LIBRARY_PATH}"
ENV CLASSPATH="${ASTRAL_PREFIX}/lib:${CLASSPATH:-}"

LABEL astral.version="${ASTRAL_COMMIT}"

COPY --from=astral_builder "${ASTRAL_PREFIX}" "${ASTRAL_PREFIX}"
COPY --from=astral_builder "${APT_REQUIREMENTS_FILE}" /build/apt/astral.txt

RUN  set -eu \
  && DEBIAN_FRONTEND=noninteractive \
  && . /build/base.sh \
  && apt-get update \
  && apt_install_from_file /build/apt/*.txt \
  && rm -rf /var/lib/apt/lists/* \
  && cat /build/apt/*.txt >> "${APT_REQUIREMENTS_FILE}"
