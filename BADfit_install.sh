#!/bin/bash
echo "Script will install BADfit and append it to your bashrc path"

echo "make sure you are in your home directory"
pwd
echo "install here[y/n] " 

read answer

if [ $answer == 'y' ]; then
   unzip BADfit.zip
   cd BADfit/bin/
   chmod +x BADfit.py 
   cd ../../
   path=`pwd`
   cd
   echo "" >> .bashrc
   echo "#added by BADfit" >> .bashrc
   echo "export PATH=\"$path/BADfit/bin:\$PATH\"" >> .bashrc
   echo "export BAD_mat=\"$path/BADfit/matlab/\"" >> .bashrc
   echo "export BAD_one=\"$path/BADfit/onetep/\"" >> .bashrc
   echo "BADfit installed"
   
else
 echo "exiting"
 
fi
    



