#!/bin/bash
function update_phi_examples {
  echo "Update Phi examples"
  cd phi/
  echo "Update .spe (single spectrum) examples"
  dataconverter SnO2_10nm.spe eln_data_phi.yaml --reader xps --nxdl NXmpes --output SnO2_10nm.spe.nxs
  echo
  echo "Update .pro (depth_profiling) examples"
  dataconverter SnO2_10nm_1.pro eln_data_phi.yaml --reader xps --nxdl NXmpes --output SnO2_10nm_1.pro.nxs
  cd ..
  echo
}

function update_scienta_examples {
  echo "Update scienta examples"
  cd scienta/
  echo "Update .txt example"
  dataconverter Cu-HHTP_*.txt eln_data.yaml --reader xps --nxdl NXmpes --output Cu-HHTP.nxs
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

update_phi_examples
update_scienta_examples
update_sle_examples