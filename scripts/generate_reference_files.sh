#!/bin/bash
function update_ref_file {
  local FOLDER=$1
  local NXDL=$2
  cd "$FOLDER"
  if [[ "$FOLDER" == "." ]]; then
    FOLDER="test"
  fi
  echo "Update $FOLDER reference file for $NXDL"
  files=$(find . -type f \( ! -name "*.log" -a ! -name "*.nxs" \))
  dataconverter ${files[@]} --reader $READER --nxdl $NXDL --ignore-undocumented --output "${FOLDER}_ref.nxs" # &> ref_output.txt
  cd ..
}

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

READER="xps"
nxdls=(
  # "NXmpes"
  "NXxps"
)

project_dir=$(dirname $(dirname $(realpath $0)))
cd $project_dir/tests/data

for folder in "${folders[@]}"; do
  for nxdl in "${nxdls[@]}"; do
    update_ref_file "$folder" "$nxdl"
  done
done