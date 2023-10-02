# to run on local with admin password
rsync -avuz --progress ~/disks/meso_S/data/amblyo_prf/derivatives/webgl/ admin@invibe.nohost.me:/var/www/my_webapp__5/www/

# to run on invibe with admin password
ssh admin@invibe.nohost.me chmod -Rfv 777 /var/www/my_webapp__5/