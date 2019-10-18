PWD=$(shell pwd)
DEBIAN_VERSION=buster-20190910-slim
TAG=snptree-v0.0.1
DOCKER_USER=darcyabjones
BASE_IMAGE=$(DOCKER_USER)/base:$(TAG)

INSTALL_BASE=/opt

define make_image_name
$(DOCKER_USER)/$(notdir $1):$(TAG)
endef

# SOFTWARE_VERSIONs

HTSLIB_TAG=1.9
BCFTOOLS_TAG=1.9
SAMTOOLS_TAG=1.9
HTSLIB_REPO=https://github.com/samtools/htslib.git
BCFTOOLS_REPO=https://github.com/samtools/bcftools.git
SAMTOOLS_REPO=https://github.com/samtools/samtools.git
HTSLIB_PREFIX_ARG=$(INSTALL_BASE)/htslib/$(HTSLIB_TAG)
BCFTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bcftools/$(BCFTOOLS_TAG)
SAMTOOLS_PREFIX_ARG=$(INSTALL_BASE)/samtools/$(SAMTOOLS_TAG)

BEDTOOLS_VERSION=2.28.0
BEDTOOLS_URL=https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz
BEDTOOLS_PREFIX_ARG=$(INSTALL_BASE)/bedtools/$(BEDTOOLS_VERSION)

SNPEFF_VERSION=4_3t
SNPEFF_URL=https://jaist.dl.sourceforge.net/project/snpeff/snpEff_v4_3t_core.zip
SNPEFF_PREFIX_ARG=$(INSTALL_BASE)/snpeff/$(SNPEFF_VERSION)

IQTREE_VERSION=1.6.12
IQTREE_URL=https://github.com/Cibiv/IQ-TREE/archive/v1.6.12.tar.gz
IQTREE_PREFIX_ARG=$(INSTALL_BASE)/iqtree/$(IQTREE_VERSION)

RAXML_COMMIT=66ad9d2
RAXML_REPO=https://github.com/amkozlov/raxml-ng.git
RAXML_PREFIX_ARG=$(INSTALL_BASE)/raxml/$(RAXML_COMMIT)

MODELTEST_VERSION=180a138
MODELTEST_REPO=https://github.com/ddarriba/modeltest.git
MODELTEST_PREFIX_ARG=$(INSTALL_BASE)/modeltest/$(MODELTEST_VERSION)

BUSCO_COMMIT=1554283ab8ee7dd5b5290f4f748234f456c36e66
BUSCO_REPO=https://gitlab.com/ezlab/busco.git
BUSCO_PREFIX_ARG=$(INSTALL_BASE)/busco/$(BUSCO_COMMIT)

AUGUSTUS_COMMIT=8b1b14a7489e4545e89c8725dc33268f6c2a9117
AUGUSTUS_REPO=https://github.com/Gaius-Augustus/Augustus.git
AUGUSTUS_PREFIX_ARG=$(INSTALL_BASE)/augustus/$(AUGUSTUS_COMMIT)

DECIPHER_VERSION=2.10.0
DECIPHER_URL=http://www2.decipher.codes/DECIPHER_2.10.0.tar.gz
DECIPHER_PREFIX_ARG=$(INSTALL_BASE)/decipher/$(DECIPHER_VERSION)

OCCULTERCUT_VERSION=1.1
OCCULTERCUT_URL=https://downloads.sourceforge.net/project/occultercut/OcculterCut_v1.1.tar.gz
OCCULTERCUT_PREFIX_ARG=$(INSTALL_BASE)/occultercut/$(OCCULTERCUT_VERSION)

ASTRAL_VERSION=5.6.3
ASTRAL_URL=https://github.com/smirarab/ASTRAL/raw/master/Astral.5.6.3.zip
ASTRAL_PREFIX_ARG=$(INSTALL_BASE)/astral/$(ASTRAL_VERSION)
ASTRAL_JAR_ARG=$(ASTRAL_PREFIX_ARG)/astral.5.6.3.jar

NWUTILS_VERSION=1.5.0
NWUTILS_URL=http://cegg.unige.ch/pub/newick-utils-1.5.0.tar.gz
NWUTILS_PREFIX_ARG=$(INSTALL_BASE)/nwutils/$(NWUTILS_VERSION)

GENERAL=bedtools python3 posix htslib snpeff busco decipher occultercut nwutils
PHYLOGENETICS=iqtree raxml astral

OPEN=$(GENERAL) $(PHYLOGENETICS)

all: docker/all singularity/all

# Docker targets
docker: docker/snptree
docker/all: docker/general docker/phylogenetics

docker/general: $(addprefix docker/,$(GENERAL))
docker/phylogenetics: $(addprefix docker/,$(PHYLOGENETICS))

# Singularity targets
singularity: singularity/snptree.sif
singularity/all: singularity/general singularity/phylogenetics

singularity/general: $(addprefix singularity/,$(addsuffix .sif,$(GENERAL)))
singularity/phylogenetics: $(addprefix singularity/,$(addsuffix .sif,$(PHYLOGENETICS)))

# Pull containers for dockerhub
dockerpull: dockerpull/snptree
dockerpull/all: dockerpull/general dockerpull/alignment dockerpull/phylogenetics
dockerpull/general: $(addprefix dockerpull/,$(GENERAL))
dockerpull/phylogenetics: $(addprefix dockerpull/,$(PHYLOGENETICS))


clean:
	docker stop $(shell docker ps -aq) || true
	docker rm $(shell docker ps -a -q) || true
	docker rmi $(shell docker images -q) || true
	singularity cache clean -a
	rm -rf -- tmp
	rm -rf -- docker

tidy:
	# Just removes untagged images to free up space
	rm -rf -- tmp
	docker image prune -f

# Generic rules
singularity/%.sif: docker/%
	@mkdir -p "$(dir $@)"
	@mkdir -p tmp tmp/cache
	SINGULARITY_TMPDIR="${PWD}/tmp" SINGULARITY_CACHEDIR="${PWD}/tmp/cache" singularity build \
	  --force \
	  --tmpdir "${PWD}/tmp" \
	  "$@" \
	  "docker-daemon://$(call make_image_name,$(notdir $(basename $@)))"

docker/%.tar.gz: docker/%
	@mkdir -p "$(dir $@)"
	docker save "$(call make_image_name,$(basename $(basename $@)))" | gzip > $@

dockerpull/%:
	docker pull "$(call make_image_name,$(notdir $@))"
	@mkdir -p "$(dir $@)"
	@touch docker/$(notdir $@)


docker/base: base.Dockerfile
	docker build \
	  --build-arg DEBIAN_VERSION=$(DEBIAN_VERSION) \
	  --build-arg MAINTAINER="$(MAINTAINER)" \
	  --build-arg PIPELINE_VERSION="$(TAG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/astral: astral.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg ASTRAL_VERSION="$(ASTRAL_VERSION)" \
	  --build-arg ASTRAL_URL="$(ASTRAL_URL)" \
	  --build-arg ASTRAL_PREFIX_ARG="$(ASTRAL_PREFIX_ARG)" \
	  --build-arg ASTRAL_JAR_ARG="$(ASTRAL_JAR_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/htslib: htslib.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_REPO="$(HTSLIB_REPO)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_REPO="$(SAMTOOLS_REPO)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --build-arg BCFTOOLS_TAG="$(BCFTOOLS_TAG)" \
	  --build-arg BCFTOOLS_REPO="$(BCFTOOLS_REPO)" \
	  --build-arg BCFTOOLS_PREFIX_ARG="$(BCFTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/python3: python3.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/posix: posix.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/bedtools: bedtools.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BEDTOOLS_VERSION="$(BEDTOOLS_VERSION)" \
	  --build-arg BEDTOOLS_URL="$(BEDTOOLS_URL)" \
	  --build-arg BEDTOOLS_PREFIX_ARG="$(BEDTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/nwutils: nwutils.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg NWUTILS_VERSION="$(NWUTILS_VERSION)" \
	  --build-arg NWUTILS_URL="$(NWUTILS_URL)" \
	  --build-arg NWUTILS_PREFIX_ARG="$(NWUTILS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/snpeff: snpeff.Dockerfile scripts/SnpSift.py docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg SNPSIFT_SCRIPT=scripts/SnpSift.py \
	  --build-arg SNPEFF_VERSION="$(SNPEFF_VERSION)" \
	  --build-arg SNPEFF_URL="$(SNPEFF_URL)" \
	  --build-arg SNPEFF_PREFIX_ARG="$(SNPEFF_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/decipher: decipher.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg DECIPHER_VERSION="$(DECIPHER_VERSION)" \
	  --build-arg DECIPHER_URL="$(DECIPHER_URL)" \
	  --build-arg DECIPHER_PREFIX_ARG="$(DECIPHER_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/iqtree: iqtree.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg IQTREE_VERSION="$(IQTREE_VERSION)" \
	  --build-arg IQTREE_URL="$(IQTREE_URL)" \
	  --build-arg IQTREE_PREFIX_ARG="$(IQTREE_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"


docker/raxml: raxml.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg RAXML_COMMIT="$(RAXML_COMMIT)" \
	  --build-arg RAXML_REPO="$(RAXML_REPO)" \
	  --build-arg RAXML_PREFIX_ARG="$(RAXML_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"


docker/modeltest: modeltest.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg MODELTEST_COMMIT="$(MODELTEST_COMMIT)" \
	  --build-arg MODELTEST_REPO="$(MODELTEST_REPO)" \
	  --build-arg MODELTEST_PREFIX_ARG="$(MODELTEST_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"


docker/augustus: augustus.Dockerfile docker/base docker/htslib
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_REPO="$(AUGUSTUS_REPO)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --build-arg HTSLIB_IMAGE="$(call make_image_name,htslib)" \
	  --build-arg HTSLIB_TAG="$(HTSLIB_TAG)" \
	  --build-arg HTSLIB_REPO="$(HTSLIB_REPO)" \
	  --build-arg SAMTOOLS_TAG="$(SAMTOOLS_TAG)" \
	  --build-arg SAMTOOLS_REPO="$(SAMTOOLS_REPO)" \
	  --build-arg HTSLIB_PREFIX_ARG="$(HTSLIB_PREFIX_ARG)" \
	  --build-arg SAMTOOLS_PREFIX_ARG="$(SAMTOOLS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@


docker/busco: busco.Dockerfile docker/base docker/augustus
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg BUSCO_COMMIT="$(BUSCO_COMMIT)" \
	  --build-arg BUSCO_REPO="$(BUSCO_REPO)" \
	  --build-arg BUSCO_PREFIX_ARG="$(BUSCO_PREFIX_ARG)" \
	  --build-arg AUGUSTUS_IMAGE="$(call make_image_name,augustus)" \
	  --build-arg AUGUSTUS_COMMIT="$(AUGUSTUS_COMMIT)" \
	  --build-arg AUGUSTUS_PREFIX_ARG="$(AUGUSTUS_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@touch $@


docker/occultercut: occultercut.Dockerfile docker/base
	docker build \
	  --build-arg IMAGE="$(BASE_IMAGE)" \
	  --build-arg OCCULTERCUT_VERSION="$(OCCULTERCUT_VERSION)" \
	  --build-arg OCCULTERCUT_URL="$(OCCULTERCUT_URL)" \
	  --build-arg OCCULTERCUT_PREFIX_ARG="$(OCCULTERCUT_PREFIX_ARG)" \
	  --file $< \
	  -t "$(call make_image_name,$@)" \
	  .
	@mkdir -p "$(dir $@)"
	@touch $@
