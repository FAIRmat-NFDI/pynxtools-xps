#!/bin/bash
READER=xps

# Function to update log
function update_log {
  local FOLDER=$1
  local NXDL=$2
  local lowercase_NXDL=$(echo "$NXDL" | tr '[:upper:]' '[:lower:]')
  log_filename="${FOLDER}_${lowercase_NXDL}_ref.log"
  echo "Generating log file at $log_filename..."
  python -c "
import os
from pynxtools.testing.nexus_conversion import get_log_file
folder = os.path.join(os.getcwd(), 'tests', 'data', '$FOLDER')
nxs_filepath = os.path.join(folder,'output.nxs')
log_filepath = os.path.join(folder,'$log_filename')
get_log_file(nxs_filepath, log_filepath, './')
"
  echo "Done!"
  echo
}

function update_log_file {
  local FOLDER=$1
  local NXDL=$2
  cd $FOLDER || exit
  echo "Update $FOLDER reference log for $NXDL"
  files=$(find . -type f \( ! -name "*.log" -a ! -name "*.nxs" \))
  dataconverter ${files[@]} --reader $READER --nxdl $NXDL --ignore-undocumented
  cd ../../.. || exit
  update_log "$FOLDER" "$NXDL"
  rm "tests/data/$FOLDER/output.nxs"
}

project_dir=$(dirname $(dirname $(realpath $0)))

folders=(
  "phi_pro"
  "phi_spe"
  "scienta_ibw"
  "scienta_txt"
  "specs_sle"
  "specs_xml"
  "specs_xy"
  "vms_irregular"
  "vms_regular"
  "vms_txt_export"
  "vms_analysis"
)

nxdls=(
  "NXmpes"
  "NXxps"
)

for folder in "${folders[@]}"; do
  for nxdl in "${nxdls[@]}"; do
    cd $project_dir/tests/data
    update_log_file "$folder" "$nxdl"
  done
done
