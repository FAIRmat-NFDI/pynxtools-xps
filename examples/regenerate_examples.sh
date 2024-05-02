#!/bin/bash
READER=xps
NXDL=NXmpes

function update_phi_examples {
  echo "Update Phi examples"
  cd phi/
  echo "Update .spe (single spectrum) example"
  dataconverter SnO2_10nm.spe eln_data_phi.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm.spe.nxs
  echo
  echo "Update .pro (depth_profiling) example"
  dataconverter SnO2_10nm_1.pro eln_data_phi.yaml --reader $READER --nxdl $NXDL --output SnO2_10nm_1.pro.nxs
  cd ..
  echo
}

function update_scienta_examples {
  echo "Update scienta examples"
  cd scienta/
  echo "Update .txt example"
  dataconverter Cu-HHTP_*.txt eln_data.yaml --reader $READER --nxdl $NXDL --output Cu-HHTP.nxs
  cd ..
  echo
}

function update_sle_examples {
  echo "Update SPECS examples"
  cd sle/
  echo "Update .sle example"
  dataconverter --params-file params.yaml
  cd ..
  echo
}

function update_vms_examples {
  echo "Update VAMAS examples"
  cd vms/
  echo "Update REGULAR file conversion example"
  dataconverter regular.vms eln_data_vms.yaml --reader $READER --nxdl $NXDL --output vms_regular_example.nxs
  echo
  echo "Update REGULAR file conversion example"
  dataconverter irregular.vms eln_data_vms.yaml --reader $READER --nxdl $NXDL --output vms_irregular_example.nxs
  cd ..
  echo
}


update_phi_examples
update_scienta_examples
update_sle_examples