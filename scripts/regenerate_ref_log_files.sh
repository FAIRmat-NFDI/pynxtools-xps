#!/bin/bash
READER=xps

#!/bin/bash

# Function to generate the filename
#!/bin/bash

# Function to generate filename
function generate_filename {
  local NXDL="$2"
  local lowercase_NXDL=$(echo "$NXDL" | tr '[:upper:]' '[:lower:]')
  local filename="$1_${lowercase_NXDL}_ref.log"
  echo "$filename"
}

# Function to update log
function update_log {
  log_filename=$(generate_filename $1 $2)
  python -c "
import os
from pynxtools.testing.nexus_conversion import get_log_file
tmp_path = os.getcwd()
get_log_file('output.nxs', '$log_filename', tmp_path)
"
}

# # Function to update phi_pro_ref
# function update_phi_pro_ref {
#   FOLDER="phi_pro"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update .pro (PHI depth_profiling) NeXus reference log for $NXDL"
#   dataconverter SnO2_10nm_1.pro eln_data_phi_pro.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_phi_spe_ref {
#   FOLDER="phi_spe"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update .spe (PHI single spectrum) NeXus reference log for $NXDL"
#   dataconverter SnO2_10nm.spe eln_data_phi_spe.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_scienta_ibw_ref {
#   FOLDER="scienta_ibw"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update scienta .ibw NeXus reference log for $NXDL"
#   dataconverter Ag_000*.ibw eln_data_scienta_ibw.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_scienta_txt_ref {
#   FOLDER="scienta_txt"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update scienta txt NeXus reference log for $NXDL"
#   dataconverter Ag_000*.txt eln_data_scienta_txt.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_specs_sle_ref {
#   FOLDER="specs_sle"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update SPECS .sle NeXus reference log for $NXDL"
#   dataconverter EX439_S718_Au.sle eln_data_specs_sle.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_specs_xml_ref {
#   FOLDER="specs_xml"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update SPECS XML NeXus reference log for $NXDL"
#   dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_specs_xml.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_specs_xy_ref {
#   FOLDER="specs_xy"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update SPECS XY NeXus reference log for $NXDL"
#   dataconverter MgFe2O4_small.xy eln_data_specs_xy.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_vms_irregular_ref {
#   FOLDER="vms_irregular"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update irregular VAMAS NeXus reference log for $NXDL"
#   dataconverter irregular.vms eln_data_vms_irregular.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_vms_regular_ref {
#   FOLDER="vms_regular"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update regular VAMAS NeXus reference log for $NXDL"
#   dataconverter regular.vms eln_data_vms_regular.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function update_vms_txt_export_ref {
#   FOLDER="vms_txt_export"
#   NXDL=$1
#   cd $FOLDER || exit
#   echo "Update Vamas text export reference log for $NXDL"
#   dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader $READER --nxdl $NXDL --ignore-undocumented
#   update_log "$FOLDER" "$NXDL"
#   rm output.nxs
#   cd ..
# }

# function cli_dataconverter {
#   # Call the CLI function with all provided files and additional flags
#   dataconverter "$@" --reader $READER --nxdl $NXDL --ignore-undocumented
# }

# cli_function() {
#   # Assuming $READER and $NXDL are defined somewhere in your script
#   READER="your_reader_value"
#   NXDL="your_nxdl_value"

#   # Accumulate all file paths into an array
#   files=("$@")

#   # Call dataconverter with all files and additional flags
#   echo "Running dataconverter with files:"
#   for file in "${files[@]}"; do
#     echo "$file"
#   done

#   # Replace with your actual dataconverter command
#   dataconverter "${files[@]}" --reader "$READER" --nxdl "$NXDL" --ignore-undocumented
# }

# # Find and filter files, then pass them to cli_function


function update_log_file {
  local FOLDER=$1
  local NXDL=$2
  cd $FOLDER || exit
  echo "Update $FOLDER reference log for $NXDL"
  files=$(find . -type f \( ! -name "*.log" -a ! -name "*.nxs" \))
  dataconverter ${files[@]} --reader $READER --nxdl $NXDL --ignore-undocumented
  update_log "$FOLDER" "$NXDL"
  rm output.nxs
  cd ..
}

project_dir=$(dirname $(dirname $(realpath $0)))
cd $project_dir/tests/data

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
)

nxdls=(
  "NXmpes"
  "NXxps"
)

for folder in "${folders[@]}"; do
  for nxdl in "${nxdls[@]}"; do
    update_log_file "$folder" "$nxdl"
  done
done
