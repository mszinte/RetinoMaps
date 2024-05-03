# -----------------------------------------------------------------------------------------
# send_data.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Send the data to webapp on server
# -----------------------------------------------------------------------------------------
# Input(s):
# none
# -----------------------------------------------------------------------------------------
# Output(s):
# sent data
# -----------------------------------------------------------------------------------------
# To run:
# On invibe.nohost.me
# 1. cd to function
# >> cd cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/webgl/
# 2. run bash command
# >> sh send_data.sh
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/prf/webgl/
# sh send_data.sh
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (mail@martinszinte.net)
# Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
# -----------------------------------------------------------------------------------------

# to run on local with admin password
rsync -avuz --progress ~/disks/meso_S/data/RetinoMaps/derivatives/webgl/ admin@invibe.nohost.me:/var/www/my_webapp/www/

# to run on invibe with admin password
ssh admin@invibe.nohost.me chmod -Rfv 777 /var/www/my_webapp/