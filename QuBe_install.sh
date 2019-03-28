#!/bin/bash
echo "Script will install QUBEKit and append it to your bashrc path"

echo "make sure you are in your home directory"
pwd
echo "install here[y/n] " 

read answer

if [ $answer == 'y' ]; then
   mkdir QuBeKit
   mv bin/ matlab/ onetep/ README.md QuBeKit
   cd QuBeKit/bin
   chmod +x QuBeKit.py 
   cd ../../
   path=`pwd`
   cd
   echo "" >> .bashrc
   echo "#added by QuBeKit" >> .bashrc
   echo "export PATH=\"$path/QuBeKit/bin:\$PATH\"" >> .bashrc
   echo "export QuBeKit=\"$path/QuBeKit/\"" >> .bashrc
   echo "QuBeKit installed"
   
else
 echo "exiting"
 
fi
    



