#!/bin/bash
READER=xps
NXDL=NXmpes

function update_cs_example {
  echo "Update coordinate system example"
  cd coordinate_system/
  dataconverter config_vms_cs_fixed.json --reader xps --nxdl NXxps --output vms-cs-fixed.nxs
  cd ..
}

function update_phi_examples {
  echo "Update Phi examples"
  cd phi/
  echo "Update .spe (single spectrum) example"
  dataconverter SnO2_10nm.spe eln_data_phi.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm.spe.nxs
  echo
  echo "Update .pro (depth_profiling) example"
  dataconverter SnO2_10nm_1.pro eln_data_phi.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm_1.pro.nxs
  cd ..
}

function update_scienta_examples {
  echo "Update scienta examples"
  cd scienta/
  cd ibw/
  echo "Update .ibw example"
  dataconverter Cu-HHTP_*.ibw eln_data_scienta_ibw.yaml --reader $READER --nxdl $NXDL --output Cu-HHTP.ibw.nxs
  cd ../txt
  echo "Update .txt example"
  dataconverter Cu-HHTP_*.txt eln_data_scienta_txt.yaml --reader $READER --nxdl $NXDL --output Cu-HHTP.txt.nxs
  cd ../..
}

function update_specs_examples {
  echo "Update SPECS examples"
  cd specs/
  cd sle/
  echo "Update .sle example"
  dataconverter EX439_S718_Au.sle eln_data_sle.yaml --reader $READER --nxdl $NXDL --output Au_25_mbar_O2_no_align.nxs
  cd ../xml
  echo "Update .xml example"
  dataconverter In-situ_PBTTT_XPS_SPECS.xml eln_data_xml.yaml --reader $READER --nxdl $NXDL --output In-situ_PBTTT.nxs
  cd ../xy
  echo "Update .xy example"
  dataconverter MgFe2O4.xy eln_data_xy.yaml --reader $READER --nxdl $NXDL --output MgFe2O4.nxs
  cd ../..
}

function update_vms_examples {
  echo "Update VAMAS examples"
  cd vms/
  echo "Update REGULAR file conversion example"
  dataconverter regular.vms eln_data_vms.yaml --reader $READER --nxdl $NXDL --output regular.vms.nxs
  echo
  echo "Update REGULAR file conversion example"
  dataconverter irregular.vms eln_data_vms.yaml --reader $READER --nxdl $NXDL --output irregular.vms.nxs
  echo "Update txt export example"
  dataconverter vms_txt_export.txt eln_data_vms_txt_export.yaml --reader $READER --nxdl $NXDL --output vms_txt_export.nxs
  cd ..
}

project_dir=$(dirname $(dirname $(realpath $0)))
cd $project_dir/examples

update_cs_example
update_phi_examples
update_scienta_examples
update_specs_examples
update_vms_examples