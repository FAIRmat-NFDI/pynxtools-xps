#!/bin/bash
READER=xps
NXDL=NXmpes


function update_phi_pro_ref {
  cd phi_pro/
  echo "Update .pro (PHI depth_profiling) NeXus reference"
  dataconverter SnO2_10nm_1.pro eln_data_phi_pro.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output phi_pro_ref.nxs  
  cd ..
}

function update_phi_spe_ref {
  cd phi_spe/
  echo "Update .spe (PHI single spectrum) NeXus reference"
  dataconverter SnO2_10nm.spe eln_data_phi_spe.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output phi_spe_ref.nxs
  cd ..
}

function update_scienta_ibw_ref {
  cd scienta_ibw/
  echo "Update scienta .ibw NeXus reference"
  dataconverter Ag_000*.ibw eln_data_scienta_ibw.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output scienta_ibw_ref.nxs
  cd ..
}

function update_scienta_txt_ref {
  cd scienta_txt/
  echo "Update scienta txt NeXus reference"
  dataconverter Ag_000*.txt eln_data_scienta_txt.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output scienta_txt_ref.nxs
  cd ..
}

function update_specs_sle_ref {
  echo "Update SPECS .sle NeXus reference"
  cd specs_sle/
  echo "Update .sle example"
  dataconverter EX439_S718_Au.sle eln_data_specs_sle.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output specs_sle_ref.nxs
  cd ..
}

function update_specs_xml_ref {
  cd specs_xml/
  echo "Update SPECS XML NeXus reference"
  dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_specs_xml.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output specs_xml_ref.nxs
  cd ..
}

function update_specs_xy_ref {
  cd specs_xy/
  echo "Update SPECS XY NeXus reference"
  dataconverter MgFe2O4_small.xy eln_data_specs_xy.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output specs_xy_ref.nxs
  cd ..
}

function update_vms_irregular_ref {
  cd vms_irregular/
  echo "Update irregular VAMAS NeXus reference"
  dataconverter irregular.vms eln_data_vms_irregular.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output vms_irregular_ref.nxs
  cd ..
}

function update_vms_regular_ref {
  cd vms_regular/
  echo "Update regular VAMAS NeXus reference"
  dataconverter regular.vms eln_data_vms_regular.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output vms_regular_ref.nxs
  cd ..
}

function update_vms_txt_export_ref {
  cd vms_txt_export/
  echo "Update Vamas text export reference"
  dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader $READER --nxdl $NXDL --ignore-undocumented --output vms_txt_export_ref.nxs
  cd ..
}

project_dir=$(dirname $(dirname $(realpath $0)))
cd $project_dir/tests/data

update_phi_pro_ref
update_phi_spe_ref
update_scienta_ibw_ref
update_scienta_txt_ref
update_specs_sle_ref
update_specs_xml_ref
update_specs_xy_ref
update_vms_irregular_ref
update_vms_regular_ref
update_vms_txt_export_ref