#!/bin/bash/

## Sets up R (3.6) and python (3.6) environments
## Written by Eduardo Maury (eduardo_maury at hms dot harvard dot edu)

USER_ID=$(whoami)
. /home/${USER_ID}/anaconda3/etc/profile.d/conda.sh
conda create -n morst_cnv python=3.6 -y
conda activate morst_cnv
conda install -y r-essentials r-base r-compquadform