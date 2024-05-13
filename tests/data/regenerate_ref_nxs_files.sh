#!/bin/bash
READER=xps
NXDL=NXmpes


function update_phi_pro_ref {
  cd phi_pro/
  echo "Update .pro (PHI depth_profiling) NeXus reference"
  dataconverter.exe SnO2_10nm_1.pro eln_data_phi_pro.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm_1.pro_ref.nxs
  cd ..
  echo
}

function update_phi_spe_ref {
  cd phi_spe/
  echo "Update .spe (PHI single spectrum) NeXus reference"
  dataconverter SnO2_10nm.spe eln_data_phi_spe.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm.spe_ref.nxs
  cd ..
  echo
}

function update_scienta_ibw_ref {
  cd scienta_txt/
  echo "Update scienta .ibw NeXus reference"
  cd scienta_txt/
  dataconverter Ag_000*.ibw eln_data.yaml --reader $READER --nxdl $NXDL --output Ag_ref.nxs
  cd ..
  echo
}

function update_scienta_txt_ref {
  cd scienta_txt/
  echo "Update scienta txt NeXus reference"
  cd scienta_txt/
  dataconverter Ag_000*.txt eln_data.yaml --reader $READER --nxdl $NXDL --output Ag_ref.nxs
  cd ..
  echo
}

function update_specs_sle_ref {
  echo "Update SPECS .sle NeXus reference"
  cd specs_sle/
  echo "Update .sle example"
  dataconverter EX439_S718_Au.sle eln_data_sle.yaml --reader $READER --nxdl $NXDL --output SPECS_SLE_ref.nxs
  cd ..
  echo
}

function update_specs_xml_ref {
  cd specs_xml/
  echo "Update SPECS XML NeXus reference"
  dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml --reader $READER --nxdl $NXDL --output SPECS_XML_ref.nxs
  cd ..
  echo
}

function update_specs_xy_ref {
  cd specs_xy/
  echo "Update SPECS XY NeXus reference"
  dataconverter MgFe2O4_small.xy eln_data.yaml --reader $READER --nxdl $NXDL --output SPECS_XY_ref.nxs
  cd ..
  echo
}

function update_vms_irregular_ref {
  cd vms_irregular/
  echo "Update irregular VAMAS NeXus reference"
  dataconverter irregular.vms eln_data_vms_irregular.yaml --reader $READER --nxdl $NXDL --output vms_irregular_ref.nxs
  cd ..
  echo
}

function update_vms_regular_ref {
  cd vms_regular/
  echo "Update regular VAMAS NeXus reference"
  dataconverter regular.vms eln_data_vms_regular.yaml --reader $READER --nxdl $NXDL --output vms_regular_ref.nxs
  cd ..
  echo
}

update_phi_pro_ref
update_phi_spe_ref
# update_scienta_ibw_ref
update_scienta_txt_ref
update_specs_sle_ref
update_specs_xml_ref
update_specs_xy_ref
update_vms_irregular_ref
update_vms_regular_ref