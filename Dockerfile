FROM ubuntu:18.04
ARG BUILD_DATE
ARG VCS_REF
LABEL maintainer="rhancock@gmail.com"
LABEL org.label-schema.name="rhancock/auto_voxel"
LABEL org.label-schema.description="Template-based VOI placement for Siemens"
LABEL org.label-schema.vcs-url="https://github.com/rhancockn/autovps"
LABEL org.label-schema.build-date=$BUILD_DATE
LABEL org.label-schema.vcs-ref=$VCS_REF

ARG DEBIAN_FRONTEND=noninteractive

ENV LANG="C.UTF-8" \
    LC_ALL="C.UTF-8"

RUN apt-get update -y && \
	apt-get install -yq --no-install-recommends \
	curl tar ca-certificates unzip \
    python3.6 python3-pip python3-setuptools && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/python3.6 /usr/bin/python

# RUN curl -fsSL -o miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh && \
#     bash miniconda.sh -b -p /usr/local/miniconda && \
#     rm miniconda.sh
# ENV PATH="/usr/local/miniconda/bin:$PATH"

# dcm2niix
RUN curl -LO https://github.com/rordenlab/dcm2niix/releases/download/v1.0.20190902/dcm2niix_lnx.zip  && \
	unzip dcm2niix_lnx.zip && \
	mv dcm2niix /usr/local/bin && \
    rm dcm2niix_lnx.zip
ENV PATH="/usr/local/bin:$PATH"

# ROBEX
RUN curl "https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1" \
        | tar -xz -C /usr/local/
ENV PATH="/usr/local/ROBEX:$PATH"

# Install FSL
COPY fsl/ /usr/local/fsl/
# Configure environment
ENV FSLDIR=/usr/local/fsl
ENV FSL_DIR="${FSLDIR}" \
    FSLOUTPUTTYPE=NIFTI \
    PATH=${FSLDIR}/bin:$PATH \
    FSLMULTIFILEQUIT=TRUE \
    POSSUMDIR=${FSLDIR} \
    LD_LIBRARY_PATH=${FSLDIR}/lib:$LD_LIBRARY_PATH \
    FSLTCLSH=/usr/bin/tclsh \
    FSLWISH=/usr/bin/wish

COPY ./ /app/
WORKDIR /app

RUN pip3 install -r requirements.txt && \
    python setup.py install && rm -r /app
COPY siemens_auto_voxel/ /siemens_auto_voxel/
ENTRYPOINT ["/siemens_auto_voxel/run.sh"]
