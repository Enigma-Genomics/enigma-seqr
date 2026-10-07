#!/bin/bash

#
# VEP init action for dataproc
#
# adapted/copied from
# https://github.com/broadinstitute/gnomad_methods/blob/main/init_scripts/vep105-init.sh
# and gs://hail-common/hailctl/dataproc/0.2.128/vep-GRCh38.sh
# 
# NB: This is code used for initializing a dataproc cluster and runs as an intialization
# action when the rest of our code is unavailable. 
#

set -euo pipefail

export PROJECT="$(gcloud config get-value project)"
export DEPLOYMENT_TYPE="$(/usr/share/google/get_metadata_value attributes/DEPLOYMENT_TYPE)"
export REFERENCE_GENOME="$(/usr/share/google/get_metadata_value attributes/REFERENCE_GENOME)"
export PIPELINE_RUNNER_APP_VERSION="$(/usr/share/google/get_metadata_value attributes/PIPELINE_RUNNER_APP_VERSION)"
export REFERENCE_DATASETS_DIR="$(/usr/share/google/get_metadata_value attributes/REFERENCE_DATASETS_DIR)"
export PIPELINE_RUNNER_BUILD_BASE="$(/usr/share/google/get_metadata_value attributes/PIPELINE_RUNNER_BUILD_BASE)"
export VEP_IMAGE_URI="$(/usr/share/google/get_metadata_value attributes/VEP_IMAGE_URI)"
export VEP_UTR_PLUGIN_SHA256="$(/usr/share/google/get_metadata_value attributes/VEP_UTR_PLUGIN_SHA256)"
[[ "$PIPELINE_RUNNER_APP_VERSION" =~ ^[0-9a-f]{40}$ ]]
[[ "$PIPELINE_RUNNER_BUILD_BASE" == gs://* ]]
if [[ "$REFERENCE_GENOME" == GRCh38 && -n "$VEP_IMAGE_URI" ]]; then
  [[ "$VEP_IMAGE_URI" =~ ^[a-zA-Z0-9./:_-]+@sha256:[0-9a-f]{64}$ ]]
  [[ "$VEP_UTR_PLUGIN_SHA256" =~ ^[0-9a-f]{64}$ ]]
fi


# Install docker
apt-get update
apt-get -y install \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg2 \
    software-properties-common \
    tabix
curl -fsSL https://download.docker.com/linux/debian/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/debian $(lsb_release -cs) stable"
apt-get update
apt-get install -y --allow-unauthenticated docker-ce

# https://github.com/hail-is/hail/issues/12936
sleep 60
sudo service docker restart

cat >/vep.c <<EOF
#include <unistd.h>
#include <stdio.h>

int
main(int argc, char *const argv[]) {
  if (setuid(geteuid()))
    perror( "setuid" );

  execv("/vep.bash", argv);
  return 0;
}
EOF
gcc -Wall -Werror -O2 /vep.c -o /vep
chmod u+s /vep

gcloud storage cp "$PIPELINE_RUNNER_BUILD_BASE/$DEPLOYMENT_TYPE/$PIPELINE_RUNNER_APP_VERSION/bin/download_vep_reference_data.bash" /download_vep_reference_data.bash
chmod +x /download_vep_reference_data.bash
./download_vep_reference_data.bash "$REFERENCE_GENOME"

gcloud storage cp "$PIPELINE_RUNNER_BUILD_BASE/$DEPLOYMENT_TYPE/$PIPELINE_RUNNER_APP_VERSION/bin/vep" /vep.bash
chmod +x /vep.bash


# Pin runtime selection for the setuid wrapper and verify the exact plugin source.
if [[ "$REFERENCE_GENOME" == GRCh38 && -n "$VEP_IMAGE_URI" ]]; then
  VEP_IMAGE_HOST=${VEP_IMAGE_URI%%/*}
  gcloud auth configure-docker "$VEP_IMAGE_HOST" --quiet
  docker pull "$VEP_IMAGE_URI"
ACTUAL_PLUGIN_SHA=$(docker run --rm --entrypoint sha256sum "$VEP_IMAGE_URI" /plugins/UTRAnnotator.pm | cut -d ' ' -f 1)
test "$ACTUAL_PLUGIN_SHA" = "$VEP_UTR_PLUGIN_SHA256"
printf '%s\n' "$VEP_IMAGE_URI" > /etc/seqr-vep-image
  chmod 0644 /etc/seqr-vep-image
fi
# Use the versioned parser contract, leaving shared reference objects unchanged.
gcloud storage cp "$PIPELINE_RUNNER_BUILD_BASE/$DEPLOYMENT_TYPE/$PIPELINE_RUNNER_APP_VERSION/vep/$REFERENCE_GENOME/vep-$REFERENCE_GENOME.json" "/var/seqr/vep-reference-data/$REFERENCE_GENOME/vep-$REFERENCE_GENOME.json"
