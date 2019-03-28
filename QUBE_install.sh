#!/bin/bash
echo "Script will install QUBEKit and append it to your bashrc path"

echo "make sure you are in your home directory"
pwd
echo "install here[y/n] " 

read answer

if [ $answer == 'y' ]; then
   mkdir QUBEKit
   mv bin/ matlab/ onetep/ README.md QUBEKit
   cd QUBEKit/bin
   chmod +x QUBEKit.py 
   cd ../../
   path=`pwd`
   cd
   echo "" >> .bashrc
   echo "#added by QUBEKit" >> .bashrc
   echo "export PATH=\"$path/QUBEKit/bin:\$PATH\"" >> .bashrc
   echo "export QUBEKit=\"$path/QUBEKit/\"" >> .bashrc
   echo "QUBEKit installed"
   
else
 echo "exiting"
 
fi
    



