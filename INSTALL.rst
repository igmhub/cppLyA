Installation
-------------------

git clone git@gitlab.in2p3.fr:jguy/baocorr.git
cd baocorr

./waf configure --prefix=..

./waf install

on lpnp204

export PKG_CONFIG_PATH=/home/guy/software/DESI/lib/pkgconfig

./waf configure --prefix=/home/guy/software/baocorr/

./waf install


Utilisation de git
-------------------
git clone git@gitlab.in2p3.fr:jguy/baocorr.git

cd baocorr

git branch toto

git checkout toto

# modifications ...

git commit -m "explication" fichier

git commit -m "explication" -a

# verif

git status

# pousser sur gitlab

git push origin toto



