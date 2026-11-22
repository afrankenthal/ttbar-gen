# Same software stack as madlad-hyper.def, as a Docker image, for building where real
# root is available (GitHub Actions) when the cluster cannot build containers.
# On the cluster:  apptainer pull madlad-hyper.sif docker://ghcr.io/<owner>/madlad-hyper:latest
#
# Versions follow the HyPER dataset record: MG5_aMC 2.9.16, Pythia 8.306, Delphes 3.5.0,
# FastJet 3.4.0, LHAPDF 6.3.0, NNPDF3.0nlo. NOT BUILD-TESTED when written.
FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive \
    LC_ALL=C \
    ROOTSYS=/opt/root \
    PATH=/opt/root/bin:/usr/local/bin:/usr/local/sbin:/usr/sbin:/usr/bin:/sbin:/bin \
    LD_LIBRARY_PATH=/usr/local/lib:/opt/root/lib:/usr/lib:/lib \
    PYTHONPATH=/opt/root/lib \
    CMAKE_PREFIX_PATH=/opt/root \
    PYTHIA8DATA=/usr/local/share/Pythia8/xmldoc

RUN apt-get update -y && apt-get install -y binutils cmake dpkg-dev g++ gcc gfortran git wget rsync unzip bc \
        ghostscript nano libssl-dev libx11-dev libxext-dev libxft-dev libxpm-dev \
        libtbb-dev libgif-dev libpcre3-dev libglu1-mesa-dev libglew-dev libftgl-dev \
        libfftw3-dev libcfitsio-dev libgsl-dev libxml2-dev liblzma-dev libxxhash-dev \
        liblz4-dev libzstd-dev zlib1g-dev python3 python3-dev python3-pip python3-six \
        python3-numpy \
    && ln -fs /usr/bin/python3 /usr/bin/python \
    && pip3 install "uproot>=5" "awkward>=2" \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ROOT (binary release for Ubuntu 20.04). If this URL has moved, pick any ROOT 6 binary
# for ubuntu20 from https://root.cern/install/all_releases/
RUN cd /opt \
    && wget -q https://root.cern/download/root_v6.26.14.Linux-ubuntu20-x86_64-gcc9.4.tar.gz \
    && tar -xzf root_v6.26.14.Linux-ubuntu20-x86_64-gcc9.4.tar.gz \
    && rm root_v6.26.14.Linux-ubuntu20-x86_64-gcc9.4.tar.gz \
    && root-config --version

# HepMC2
RUN mkdir -p /opt/build && cd /opt/build \
    && wget -q https://hepmc.web.cern.ch/hepmc/releases/HepMC-2.06.11.tar.gz \
    && tar -zxf HepMC-2.06.11.tar.gz && mkdir HepMC_build && cd HepMC_build \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local -Dmomentum=GEV -Dlength=MM ../HepMC-2.06.11 \
    && make -j4 && make install && rm -rf /opt/build

# LHAPDF 6.3.0 and NNPDF3.0 NLO (LHAPDF id 260000)
RUN mkdir -p /opt/build && cd /opt/build \
    && wget -q "https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz" -O LHAPDF-6.3.0.tar.gz \
    && tar -zxf LHAPDF-6.3.0.tar.gz && cd LHAPDF-6.3.0 \
    && ./configure --prefix=/usr/local && make -j4 && make install && rm -rf /opt/build \
    && cd /usr/local/share/LHAPDF \
    && wget -q http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF30_nlo_as_0118.tar.gz \
    && tar -zxf NNPDF30_nlo_as_0118.tar.gz && rm NNPDF30_nlo_as_0118.tar.gz

# Pythia 8.306 (MadLAD fetches .tar; the site also serves .tgz)
RUN mkdir -p /opt/build && cd /opt/build \
    && ( wget -q https://pythia.org/releases/pythia83/pythia8306.tgz \
         || wget -q https://pythia.org/download/pythia83/pythia8306.tgz ) \
    && tar -xf pythia8306.t* && cd pythia8306 \
    && ./configure --prefix=/usr/local --with-hepmc2=/usr/local \
        --with-hepmc2-include=/usr/local/include --with-lhapdf6=/usr/local \
        --with-lhapdf6-plugin=LHAPDF6.h --with-gzip=/usr \
    && make -j4 && make install && rm -rf /opt/build

# FastJet 3.4.0
RUN mkdir -p /opt/build && cd /opt/build \
    && wget -q http://fastjet.fr/repo/fastjet-3.4.0.tar.gz \
    && tar -zxf fastjet-3.4.0.tar.gz && cd fastjet-3.4.0 \
    && ./configure --prefix=/usr/local && make -j4 && make install && rm -rf /opt/build

# Delphes 3.5.0
RUN mkdir -p /opt/build && cd /opt/build \
    && wget -q http://cp3.irmp.ucl.ac.be/downloads/Delphes-3.5.0.tar.gz \
    && tar -zxf Delphes-3.5.0.tar.gz && mkdir Delphes_build && cd Delphes_build \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local ../Delphes-3.5.0 \
    && make -j4 install && rm -rf /opt/build

# MG5_aMC 2.9.16, configured as MadLAD does. The trial NLO process at the end makes MG5
# install its loop libraries now; the image is read-only at run time.
RUN mkdir -p /app/MG5_aMC /opt/build && cd /opt/build \
    && wget -q https://launchpad.net/mg5amcnlo/lts/2.9.x/+download/MG5_aMC_v2.9.16.tar.gz \
    && tar -zxf MG5_aMC_v2.9.16.tar.gz -C /app/MG5_aMC --strip-components 1 \
    && cd /app/MG5_aMC \
    && ln -fs /app/MG5_aMC/bin/mg5_aMC /usr/local/bin/mg5 \
    && sed -i 's/# pythia8_path = .\/HEPTools\/pythia8/pythia8_path = \/usr\/local/' ./input/mg5_configuration.txt \
    && sed -i 's/# delphes_path = .\/Delphes/delphes_path = \/usr\/local/' ./input/mg5_configuration.txt \
    && sed -i 's/# lhapdf_py3 = lhapdf-config/lhapdf_py3 = \/usr\/local\/bin\/lhapdf-config/' ./input/mg5_configuration.txt \
    && sed -i 's/# fastjet = fastjet-config/fastjet = \/usr\/local\/bin\/fastjet-config/' ./input/mg5_configuration.txt \
    && sed -i 's/# auto_update = 7/auto_update = 0/' ./input/mg5_configuration.txt \
    && cd /opt/build \
    && printf 'import model loop_sm-no_b_mass\ninstall mg5amc_py8_interface\ngenerate p p > t t~ [QCD]\noutput mg5_test_run\n' > mg5_exec_card \
    && /usr/local/bin/mg5 mg5_exec_card \
    && rm -rf /opt/build

# sanity checks: fail the build here rather than in the first job
RUN which DelphesHepMC2 && pythia8-config --version && lhapdf-config --version && fastjet-config --version \
    && ls /usr/local/share/LHAPDF/NNPDF30_nlo_as_0118 >/dev/null \
    && python3 -c "import uproot, awkward; print('uproot', uproot.__version__)"
