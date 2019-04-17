FROM multiarch/ubuntu-debootstrap:arm64-bionic

WORKDIR /build
COPY . .

RUN uname -a
RUN apt-get update -qq && \
  apt-get install -yq --no-install-suggests --no-install-recommends \
  build-essential \
  ca-certificates \
  gcc \
  git \
  g++ \
  make \
  software-properties-common \
  zlib1g-dev
RUN add-apt-repository -y universe && \
  apt-get install -yq \
  libtbb-dev

CMD bash
