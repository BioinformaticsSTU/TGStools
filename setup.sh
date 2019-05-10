#!/bin/sh
#python 
tar zvxf PLEK.1.2.tar.gz 
cd PLEK.1.2
python PLEK_setup.py
cd ../
unzip CNCI-master.zip
cd CNCI-master
chmod ugo+x  twoBitToFa
unzip libsvm-3.0.zip
cd libsvm-3.0
make
cd ../../
conda install -y numpy='1.13.3'
conda install -y pandas='0.24.1'
conda install -y matplotlib='3.0.3'
conda install -y gseapy='0.9.9'
conda install -y scipy='1.2.1'



