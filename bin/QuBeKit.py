#!/usr/bin/env python
#==================================================================================================
#|      _____                ____               __  __             __                             |
#|     /\  __`\             /\  _`\            /\ \/\ \     __    /\ \__                          |
#|     \ \ \/\ \    __  __  \ \ \L\ \     __   \ \ \/'/'   /\_\   \ \ ,_\                         |
#|      \ \ \ \ \  /\ \/\ \  \ \  _ <'  /'__`\  \ \ , <    \/\ \   \ \ \/                         |
#|       \ \ \\'\\ \ \ \_\ \  \ \ \L\ \/\  __/   \ \ \\`\   \ \ \   \ \ \_                        |
#|        \ \___\_\ \ \____/   \ \____/\ \____\   \ \_\ \_\  \ \_\   \ \__\                       | 
#|         \/__//_/  \/___/     \/___/  \/____/    \/_/\/_/   \/_/    \/__/                       |
#==================================================================================================
#|                                 Quantum Bespoke-kit                                            |
#==================================================================================================
#Utility for the derivation of specific ligand parameters

import sys
import numpy as np
import argparse
import os
from string import digits
import math
from time import gmtime, strftime
import fileinput
import collections
import re
from decimal import Decimal
import matplotlib.pyplot as plt
remove_digits = str.maketrans('', '', digits)
from colorama import Fore
from colorama import Style
try:
    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk.unit import *
except:
    print('OpenMM not found single point comparison not available')


###################################################################################################
#Parameters to be editted
###################################################################################################
theory='wB97XD/6-311++G(d,p)' #g09 theory to use in freq and dihedral scans recomended wB97XD/6-311++G(d,p)
vib_scaling=0.957  #the associated scalling to the theory
processors=2 #the amount of processors to be used in the g09 scripts this affects the bonds and dihedral scans
memory=2 #the amount of memory in GB to be specificed in the g09 scripts
dihstart=0 #starting angle of dihedral scan
increment=15 #angle increse increment
numscan=25 # the number of optimisations around the dihedral angle
T_wieght=999999 #the wieghting temperature that can be changed to better fit complicated surfaces
improper_list=[0, 160, 161, 162, 165, 205, 221, 277] #list of improper torsions not to refit includes 165 a common ring torsion
new_dihnum=501      #the parameter number for the new dihedral to be fit 
Q_file='results.dat' #if the results are collected with QuBeKit this is always true
tor_limit=20        #torsion Vn limit to speed up fitting 
divisonarray= [2, 1, 0.5, 0.1, 0.01, 0.001] #parameter divison array 
div_index =0
XML_GMX_list=[0, 160, 161, 162, 221, 277] #this list should not need to be changed and seperates proper torsions and out of plane improper torsions for the xml and GMX files
###################################################################################################
parser = argparse.ArgumentParser(prog='QuBeKit.py', formatter_class=argparse.RawDescriptionHelpFormatter,
description="""
Utility for the derivation of specific ligand parameters
Requires BOSS make sure BOSSdir is set in bashrc
Example input to write the bond and angles input file for g09 from a zmat file
python QuBeKit.py -f bonds -t write -z toluene.z -c 0 -m 1
File names NB/NBV non-bonded with/out virtual sites, BA bonds and angles fited, D dihedrals fitted
Final file name QuBe will be wrote when xml and GMX files are made
""")
parser.add_argument('-f', "--function", help='Enter the function you wish to use bonds (covers bonds and angles terms), dihedrals and charges etc')
parser.add_argument('-t', "--type", help='Enter the function type you want this can be write , fit or analyse in the case of dihedrals (xyz charge input is wrote when bonds are fit) when writing and fitting dihedrals you will be promted for the zmat index scanning')
parser.add_argument('-X', "--XML", help='Option for making a XML file  and GMX gro and itp files if needed options yes default is no')
parser.add_argument('-z', '--zmat', help='The name of the zmat with .z')
parser.add_argument('-p', '--PDB', help='The name of the pdb file with .pdb')
parser.add_argument('-c', '--charge', help='The charge of the molecule nedded for the g09 input files, defulat = 0')
parser.add_argument('-m', '--multiplicity', help='The multiplicity of the molecule nedded for the g09 input files, defult = 1')
parser.add_argument('-s', '--submission', help='Write a sbmission script for the function called default no ')
parser.add_argument('-v', '--Vn', help='The amount of Vn coefficients to fit in the dihedral fitting, default = 4 ')
parser.add_argument('-l', '--penalty', help='The penalty used in torsion fitting between new parameters and opls reference, default = 0')
parser.add_argument('-d', '--dihedral', help='Enter the dihedral number to be fit, default will look at what SCAN folder it is running in')
parser.add_argument('-FR', '--frequency', help='Option to perform a QM MM frequency comparison, options yes default no')
parser.add_argument('-SP', '--singlepoint', help='Option to perform a single point energy calculation in openMM to make sure the eneries match (OPLS combination rule is used)')
parser.add_argument('-r', '--replace', help='Option to replace any valid dihedral terms in a molecule with QuBeKit previously optimised values. These dihedrals will be ignored in subsequent optimizations')
args = parser.parse_args()
#Check that BOSSdir is set
if 'BOSSdir' in os.environ:
     print('BOSSdir set!')
else: 
    print('$BOSSdir missing make sure BOSSdir is set in bashrc to use all features')
if not args.zmat and not args.PDB and args.type!='analyse':
  sys.exit('Zmat or PDB missing please enter')
if args.zmat:    
   molecule_name=args.zmat[:-2] #remove .z tag
if args.PDB:
   molecule_name=args.PDB[:-4] #remove .pdb tag
if args.charge:
   charge=int(args.charge) 
else:
   charge=0
if args.multiplicity:
   multiplicity=int(args.multiplicity)
else:
   multiplicity=1
if args.Vn:
   N_parameters=int(args.Vn)
else:
   N_parameters=4
if args.penalty:
   Lambda=float(args.penalty)
else:
   Lambda=0 
if args.function=='dihedrals':
  if args.type=='fit' or args.type=='check': #If dihedrals fit is called find the dihedral number to fit    
   if args.dihedral:
      dihnum=int(args.dihedral) #Checking command line first
   else:
      where=os.getcwd()
      if "FITTING" in where and "SCAN" in where: #Then looking in scan folders
          where=where[:-8]
          where=where[-2:]
          where=re.sub("\D","",where)
          dihnum=int(where)
      elif "SCAN" in where:
          where=where[-2:]
          where=re.sub("\D","",where)
          dihnum=int(where)   
###################################################################################################
#function list
###################################################################################################
def settings():
    print('''Gaussian 09 settings used: 
memory = %sGB  
processors = %s 
charge = %s 
multiplicity = %s 
theory = %s 
qm scalling = %s'''%(memory, processors, charge, multiplicity, theory, vib_scaling))
###################################################################################################
def BONDS_write():  #Write the g09 style optimisation and freq QM script from the PDB file 
    out=open('zmat.com','w+')
    out.write("""%%Mem=%sGB
%%NProcShared=%s
%%Chk=lig
# %s SCF=XQC Opt=tight freq
    
ligand
    
%s %s
"""%(memory, processors, theory, charge, multiplicity))
    pdb=open('%s.pdb'%(molecule_name),'r').read()
    Atoms=pdb.count('HETATM') #Pdb style using HETATM
    PDB=open('%s.pdb'%(molecule_name))
    lines=PDB.readlines()
    for line in lines:
      if 'HETATM' in line:
       tag= line.split()[2][:2].translate(remove_digits)
       out.write('%-2s      %9.6f    %9.6f    %9.6f\n'%(tag, float(line.split()[6]), float(line.split()[7]), float(line.split()[8])))
    out.write('\n')
    out.write('\n')
    out.write('\n')
    out.write('\n')
    out.close()
    print('%s Atoms found and zmat.com wrote'%(Atoms))
###################################################################################################
def DIHEDRALS_write(): #Write the g09 style constrained optimisation QM script from the z-mat to keep the numbering consistent
    zmat=open('%s.z'%(molecule_name),'r')
    lines=zmat.readlines()
    i=dihstart
    while i < increment*numscan:
      os.system('mkdir run_%s'%(i))
      os.chdir('run_%s'%(i))
      out=open('zmat.com','w+')
      out.write("""%%Mem=%sGB
%%NProcShared=%s
%%Chk=lig
# %s SCF=XQC Opt=(ModRedundant) 
    
ligand
    
%s %s
"""%(memory, processors, theory, charge, multiplicity))
      for (n,line) in enumerate(lines):
        if 0 < n < G1:
         if int(line.split()[0]) == int(dihnum):
           fix1=int(dihnum)-2
           fix2=int(line.split()[4])-2
           fix3=int(line.split()[6])-2
           fix4=int(line.split()[8])-2
           out.write('%-2s%7i%12.6f%5i%12.6f%5i%12.6f\n'%(str(line.split()[1])[:2].translate(remove_digits), int(line.split()[4]), float(line.split()[5]), int(line.split()[6]), float(line.split()[7]), int(line.split()[8]),i ))
         elif n == 1:
           out.write('X\n')
         elif n== 2:
           out.write('X %7i%12.6f\n'%(int(line.split()[4]) , float(line.split()[5])))
         elif n == 3:
           out.write('%-2s%7i%12.6f%5i%12.6f\n'%(str(line.split()[1])[:2].translate(remove_digits), int(line.split()[4]), float(line.split()[5]), int(line.split()[6]), float(line.split()[7])))
         else:
           out.write('%-2s%7i%12.6f%5i%12.6f%5i%12.6f\n'%(str(line.split()[1])[:2].translate(remove_digits), int(line.split()[4]), float(line.split()[5]), int(line.split()[6]), float(line.split()[7]), int(line.split()[8]), float(line.split()[9]) ))
      out.write('\n')
      out.write('%-2i %2i %2i %2i  F'%(int(fix1), int(fix2), int(fix3), int(fix4)))
      out.write('\n')
      out.write('\n')
      out.write('\n')
      out.write('\n')
      if args.submission:
         DIHEDRALS_sub() 
         os.system('sbatch dihedrals_script.sub')
      os.chdir('../')    
      i=i+increment  
###################################################################################################
def DIHEDRALS_find(): #Find all variable dihedrals in the zmat and print to screen to help the user chosse 
    G1, V1, A1, D1, Q1 = dih_tag(molecule_name)
    with open("%s.z"%(molecule_name),"r") as f:
     lines=f.readlines()
     var_dihs=np.zeros((A1-(V1+1),2))
     i=0
     for  line in lines[V1+1:A1]: #Make an array of dihedral number and OPLS type
            var_dihs[i,0]=line.split()[0]
            var_dihs[i,1]=line.split()[1]
            i=i+1
    out=open("variable_dihedrals.dat","w+") #write all dihedrals to a dat file 
    i=0 
    with open("%s.z"%(molecule_name),"r") as f:
      lines=f.readlines()
      for  line in lines[int(var_dihs[0,0]):int(var_dihs[len(var_dihs)-1,0])+1]:
        out.write("%4s %4s %4s %4s %4i\n" %(line.split()[0],line.split()[4],line.split()[6],line.split()[8],int(var_dihs[i,1])))
        i=i+1  
    out.close()
    dihedrals=open('variable_dihedrals.dat','r')
    lines=dihedrals.readlines()
    var_no=[]
    var_types=[]
    zmat_numbers=[]
    same_atoms=[]
    for line in lines:  #Check to see if the torsion is in the list of impropers we dont want to refit 
        if int(line.split()[4]) not in improper_list and int(line.split()[4]) not in var_types:
           var_types.append(int(line.split()[4]))
        if int(line.split()[4]) not in improper_list: 
          var_no.append(int(line.split()[0]))
    print('%s variable dihedrals found: %s'%(len(var_no), var_no)) #Print the dihedrals found 
    if len(var_types) == 0:
       print('No dihedrals found to scan')
       scan=input('Would you still like to write scan files?:\n> ')
       if scan == 'yes' or scan == 'y':
          for i in range(0,len(var_types)):
              for line in lines:
                  if var_types[i] == int(line.split()[4]):
                     rot_a=int(line.split()[1])
                     rot_b=int(line.split()[2])
                     for j in range(i+1,len(var_types)):
                         for line in lines:
                             if int(line.split()[4]) == var_types[j]:
                                rot_1=int(line.split()[1])
                                rot_2=int(line.split()[2])
                                if rot_a == rot_1 and rot_b == rot_2:
                                   if var_types[j] not in same_atoms:
                                      same_atoms.append(var_types[j])
                                elif rot_a == rot_2 and rot_b == rot_1:
                                     if var_types[j] not in same_atoms:
                                        same_atoms.append(var_types[j])
       else:
             sys.exit('No scan files wrote')
    print('Minimum amount of scans needed to cover all torsions = %s'%(len(var_types))) #Print the amount of OPLS dihedral types found in the zmat this is the minimum amount of scans
    print('The zmat index numbers for these torsions are:') #Print the zmat numbers of these variable dihedrals 
    if len(same_atoms)!=0:                                 
       for i in range(0,len(same_atoms)):
           var_types.remove(same_atoms[i])
    for j in range(0,len(var_types)):
        for line in lines:
            if int(line.split()[4]) == var_types[j]:
               zmat_numbers.append(int(line.split()[0])) #Finding the zmat number based on OPLS dihedral type
               break
    print(zmat_numbers)
    os.system('rm variable_dihedrals.dat')
    return zmat_numbers
###################################################################################################
def DIHEDRALS_analyse(): #Collect the results of the QM dihedral optimisations and write them to results.dat files can handle multipule SCAN folder at once 
    out=open('results.dat','w+')
    i=dihstart
    errors=0
    #grid=np.linspace(dihstart, (increment*numscan)-increment, num=numscan)
    result=[]
    SCF=0
    for x in range(numscan):
        result.append('+')
    while i < increment*numscan:
          if "run_%s"%(i) not in os.listdir("."):
            sys.exit("run_%s folder missing QM scan should be done first before fitting"%(i))   
          os.chdir('run_%s'%(i))
          if 'zmat.log' not in os.listdir("."): #If scan result missing add X
               result[int(i/increment)]='X'  
          else: 
             log=open('zmat.log','r').read()
             error=log.count('Error ')
             if error > 0:
                errors=errors+1
                result[int(i/increment)]='X'  #If error in results of QM scan add X
             else:
               #print('Folder run_%s'%(i))
               log=open('zmat.log','r')
               lines=log.readlines()
               for (n,line) in enumerate(lines):
                    if 'SCF Done:' in line: #If file is fine find the last SCF Done line
                        SCF=n
               if SCF!=0:
                  
                 for (n,line) in enumerate(lines): 
                    if n == SCF:
                       out.write("%3i %12.6f\n"%(i, float(line.split()[4])))
                    if "Normal termination" in line:
                        result[int(i/increment)]='O'  #Add the energy to the results file and add O for okay
          os.chdir('../')  
          i=i+increment
    result=np.reshape(result, (int(len(result)),1))     
    
    return errors, result
###################################################################################################
def DIHEDRALS_fit(): #Main dihedral fitting function, calls sub functions to run the fitting
    if "results.dat" not in os.listdir("."):
        print('results file missing')
    print(r"""
==================================================================================================    
==================================================================================================
|      _____                ____               __  __             __                             |
|     /\  __`\             /\  _`\            /\ \/\ \     __    /\ \__                          |
|     \ \ \/\ \    __  __  \ \ \L\ \     __   \ \ \/'/'   /\_\   \ \ ,_\                         |
|      \ \ \ \ \  /\ \/\ \  \ \  _ <'  /'__`\  \ \ , <    \/\ \   \ \ \/                         |
|       \ \ \\'\\ \ \ \_\ \  \ \ \L\ \/\  __/   \ \ \\`\   \ \ \   \ \ \_                        |
|        \ \___\_\ \ \____/   \ \____/\ \____\   \ \_\ \_\  \ \_\   \ \__\                       | 
|         \/__//_/  \/___/     \/___/  \/____/    \/_/\/_/   \/_/    \/__/                       |
==================================================================================================
|                             Quantum Bespoke Force Field kit                                    |
==================================================================================================""")
    print('''|                         Fitting enviroment specific torsions                                   |
==================================================================================================''')
    tor_list, N_torsions, param_types, tor_string = find_torsions_types(molecule_name,dihnum) #Find the torsions to be fittied based on if they pass through the rotated bond, compare opls atoms types and renumber OPLS torsions if atom type dont match but torsion type does.
    starting_zmat(molecule_name) #
    cmd_prep(dihnum)
    starting_error( Q_file, Lambda, N_torsions)
    plot()
    os.system("cp Plot Starting_plot")
    os.system("rm Plot ")
    clean()
    os.system("rm ligandqm dihzmat")
    dih_prep(molecule_name) 
    opt(Q_file, Lambda, molecule_name,div_index,divisonarray,tor_limit,N_torsions, tor_string)
    new_zmat(molecule_name)
    #ref_find(N_torsions) 
    os.system(' rm dihedral_list.dat Enviroment_dihedral_list.dat ligandqm dihzmat log')
    os.system(' cp oplsaa.par ../../QuBe_PARAMS/QuBe.par')   
    print('''==================================================================================================''')         
###################################################################################################
def clean(): #Clean the outputs of each MM scan 
    os.system('cp oplsstore oplsaa.par')
    if 'errorout' in os.listdir("."):
     os.system('rm errorout')
    if 'tempout' in os.listdir("."):
     os.system('rm tempout')
    if 'ligandmm' in os.listdir("."):
     os.system('rm ligandmm')
    if 'compare' in os.listdir("."):
     os.system('rm compare')
    if 'dih.CSV' in os.listdir("."):
     os.system('rm dih.CSV')
    return None
###################################################################################################     
def write_to_par(torsionparams, N_torsions, tor_string): #Write the updated parameters to the oplsaa.par file between scans
     #print(tor_string)
     #print(N_torsions)
     out=open('oplsaa.par', 'a')
     for i in range(N_torsions):
         no=new_dihnum+i
         out.write('%s  % 4.3f    % 4.3f    % 4.3f    % 4.3f      %s-%s-%s-%s     **NEWQuBe**\n' %(no, torsionparams[i,0], torsionparams[i,1], torsionparams[i,2], torsionparams[i,3], tor_string[i,0], tor_string[i,1], tor_string[i,2], tor_string[i,3])) 
     out.close()
     return None
###################################################################################################    
def QM_prep(Q_file): #Prep the QM results file 
    os.system(" awk '{print $2}' %s >> ligandqm"%(Q_file))
    return None
###################################################################################################
def run_boss(): #Run the dihedral driver script
    os.system('csh dihcmd > log ')
    os.system("awk '{print $4}' dih.CSV |grep -F .| sed 's/,/0/g' >> ligandmm")
    os.system('rm dihout dih.pdb dihsum')
    return None
###################################################################################################    
def cmd_prep(dihnum): #Prep the BOSS command file to read our new custom param files if needed
    f=open('dihcmd','r')
    lines=f.readlines()
    f.close()
    out=open('dihcmd','w')
    for line in lines:
     if "set atomnumber" in line:
        out.write(line[0:18]+"%s\""%(str(dihnum).zfill(6))+"            \n")
     else:
        out.write(line)
###################################################################################################
def errorcal(Lambda,dih_ref,torsionparams, T_wieght):
    K_b=0.001987
    np.set_printoptions(suppress=True, precision=12)
    mm=open('ligandmm','r')
    qm=open('ligandqm','r') #Does not change between scans maybe keep in memory?
    linesqm=( line for line in qm)
    qm_table=np.loadtxt(linesqm, usecols=(0))
    qm_min=min(qm_table)
    linesmm=(line for line in mm)
    mm_table=np.loadtxt(linesmm, usecols=(0))
    qm_table=(qm_table-qm_min) * 627.509
    table=np.array([qm_table, mm_table])
    table = table.T
    np.savetxt('compare', table,  fmt='% 8.6f       % 8.6f')
    ERRA=(mm_table-qm_table)**2
    ERRA1=ERRA * np.exp(-qm_table/(K_b*T_wieght))
    ERRS=math.sqrt(np.sum(ERRA)/(numscan))
    ERRS1=math.sqrt(np.sum(ERRA1)/(numscan))
    penalty= Lambda * np.sum(np.absolute(dih_ref-torsionparams))
    sumerror=(ERRS1 + penalty)
    return sumerror, penalty
###################################################################################################
def starting_error( Q_file, Lambda, N_torsions):
    dih_ref = ref_find(N_torsions)
    #print(N_torsions)
    torsionparams = np.zeros((N_torsions,4))
    print('''==================================================================================================
| Finding starting error ...                                                                     |
==================================================================================================''')
    out=open('starting_error.txt','w+')
    out.write("Starting Error value OPLSAA\n")
    clean()
    QM_prep(Q_file)
    run_boss()
    sumerror, penalty =errorcal(0,dih_ref,torsionparams, T_wieght)
    out.write("Error= % 2.3f\n" %(sumerror))
    out.write("Bias= % 2.3f\n" %(penalty))
    print('''| Starting error = %2.3f                                                                         |''' %(sumerror))
    print('''| Bias = %2.3f                                                                                   |''' %(penalty))
    print('''| Done!                                                                                          |
==================================================================================================''')
    return None
###################################################################################################      
def plot():
    comp=open('compare','r')
    linescomp=( line for line in comp)
    comp_table=np.loadtxt(linescomp, usecols=(0,1))
    plot=open('Plot','w+')
    plot.write('Torsion energy comparison\n')
    plot.write('Angle  QM energy   MM energy\n')
    angles=np.linspace(0, 360, 25)
    i=0
    while i < 25:
     plot.write("% 3i   % 3.6f   % 3.6f\n"%(angles[i],comp_table[i,0],comp_table[i,1]))
     i=i+1
    return None
###################################################################################################
def py_plot():
    with open('Starting_plot','r') as f:
         data = f.read()
         
    data = data.split('\n')
    
    angle=[float(row.split('    ')[0]) for row in data[2:-1]]
    QM=[float(row.split('    ')[1]) for row in data[2:-1]]
    MM=[float(row.split('    ')[2]) for row in data[2:-1]]
    
    with open('Plot','r') as p:
         data_p = p.read()
         
    data_p = data_p.split('\n')
    fitted=[float(row.split('    ')[2]) for row in data_p[2:-1]]
    QM_data, = plt.plot(angle, QM, label='QM data', marker='o', linestyle='None', color='blue')
    BOSS_data = plt.plot(angle, MM, label='Starting parameters', linestyle='--', color='red')
    Fitted_data = plt.plot(angle, fitted, label='Final parameters', color='black')
    plt.title("Relative energy as a function of torsion angle %s"%(dihnum))
    plt.xlabel("Torsion angle$^{\circ}$")
    plt.ylabel("Relative energy Kcal/mol")
    plt.legend(loc=1)
    #fig=plt.plot(angle,QM, 'bo', angle, MM, 'r--', angle, fitted, 'k')
    plt.savefig('Scan.pdf')
###################################################################################################   
def opt(Q_file, Lambda, molecule_name, div_index,divisonarray,tor_limit,N_torsions, tor_string):
    dih_ref = ref_find(N_torsions)
    if args.penalty:
       torsionparams=np.array(dih_ref)
       origvalue = np.array(dih_ref)
    else:
         torsionparams = np.zeros((N_torsions,4))
         origvalue = np.zeros((N_torsions,4))
    if Lambda > 0:
       print('''| Starting optimization using a penalty of lambda = %3.2f and torsion coefficient limit = %2i      |'''%(float(Lambda), int(tor_limit)))
    else:   
        print('''| Starting optimisation using no penalty and torsion coefficicent limit = %2i                     |'''%(tor_limit))
    results=open('resultsout','w+')
    results.write("Started % s\n"%(strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    clean()
    QM_prep(Q_file)
    print('''| RUNNING ....                                                                                   |''')    
    n=0
    while n < 1000:
     un_changed=0
     for i in range(0,N_torsions):
      for j in range(0,N_parameters): 
       clean()  
       write_to_par(torsionparams, N_torsions, tor_string)  #write the new params to par
       run_boss()                               #do first run
       sumerror, penalty =errorcal(Lambda,dih_ref,torsionparams, T_wieght) #find error using torsionparams input
       penalty1=penalty   #set the penaltys so they can be printed
       sumerror1=sumerror #set first error as sumerror1
       torsionparams[i,j] = origvalue[i,j] + divisonarray[div_index] #increase the parameter
       clean()
       write_to_par(torsionparams, N_torsions, tor_string) #write new params to par
       run_boss()
       sumerror, penalty =errorcal(Lambda,dih_ref,torsionparams, T_wieght) #find new error
       if torsionparams[i,j] >= tor_limit:
          sumerror2= sumerror + 1
          penalty2 = penalty +1
       else:
        sumerror2=sumerror #set the second set of errors
        penalty2=penalty
       torsionparams[i,j] = origvalue[i,j] - divisonarray[div_index] #decrease the parameter
       clean()
       write_to_par(torsionparams, N_torsions, tor_string) #write new params to par
       run_boss()
       sumerror, penalty =errorcal(Lambda,dih_ref,torsionparams, T_wieght) 
       sumerror3=sumerror #set the thrid set of errors
       penalty3=penalty
       tot_error=[sumerror1, sumerror2, sumerror3]
       minerror=min(tot_error)
       #print("error = % 3.3f"%(minerror))
       if minerror == sumerror1:
        torsionparams[i,j] = origvalue[i,j]
        error_pen=penalty1
        un_changed = un_changed +1
       elif minerror == sumerror2:
        origvalue[i,j] = origvalue[i,j] + divisonarray[div_index]
        torsionparams[i,j] = origvalue[i,j]
        error_pen=penalty2
       else:
        origvalue[i,j] = origvalue[i,j] - divisonarray[div_index]
        torsionparams[i,j] = origvalue[i,j]
        error_pen=penalty3
       


       results.write('Parameters % s\n'%(torsionparams))
       results.write('Error= % 3.3f\n'%(minerror))
       results.write('Penalty= % 3.3f\n'%(error_pen))
       print('''| Curent error= %3.3f Penalty= %3.3f Parameter Increment= %4.3f                                  |'''%(minerror, error_pen, float(divisonarray[div_index])), end='\r', flush=True)
     n=n+1
     if un_changed == (N_torsions * N_parameters):
        #print(un_changed)
        div_index = div_index +1
        if div_index >= len(divisonarray):
          break
    clean()
    write_to_par(torsionparams, N_torsions, tor_string)
    run_boss()
    errorcal(Lambda,dih_ref,torsionparams, T_wieght)
    plot()
    py_plot()
    results.write("Done % s\n"%(strftime("%Y-%m-%d %H:%M:%S", gmtime())))
    print('''| Done!                                                                                          |''')
    print('''| Final Error= % 4.3f                                                                            |'''%(minerror))
    print('''| Penalty= % 4.3f                                                                                |'''%(error_pen))
    print('''| Final Parameters:                                                                              |''')
    print('''|  V1      V2      V3      V4                                                                    |''')
    for i in range(0,len(torsionparams)):
      print('''| %6.3f  %6.3f  %6.3f  %6.3f                                                                 |'''%(torsionparams[i,0], torsionparams[i,1], torsionparams[i,2], torsionparams[i,3]))        
    return None
###################################################################################################   
def starting_zmat(molecule_name):
    G1, V1, A1, D1, Q1 =dih_tag(molecule_name)
    zmat=open("%s.z"%(molecule_name),"r")
    lines=zmat.readlines()
    line=lines[dihnum].split()
    C1=line[0]   #find the atom numbers of the atoms in the dihedral
    C2=line[4]
    C3=line[6]
    C4=line[8]
    var=lines[V1+1:A1]
    for line in var:
      if line.startswith('%4i'%(int(dihnum))):
       d_type=line.split()
       dr_type=d_type[1]
       c_type=d_type[3]      
    f=open("%s.z"%(molecule_name),"r")
    lines=f.readlines()
    f.close()
    out=open("dihzmat","w+")
    for (n, line) in enumerate(lines):
        if "%4i%4i%4i"%(int(dihnum),int(dr_type),int(dr_type)) not in line:
            if n == int(dihnum):
               out.write(line[:56]+'%11.6f'%(0)+line[67:])
            elif n == D1:
               out.write("%4i%4i%4i%4i%4i%4i\n                    Domain Definitions follow     (4I4)\n"%(int(C1),int(C2),int(C3),int(C4),int(dr_type),int(dr_type)))
            else:
                 out.write(line)
    return None
###################################################################################################    
def ref_find(N_torsions):    
    #print(df_type)
    dih_ref = np.zeros((N_torsions,4))
    par=open("oplsaa.par","r")
    lines = par.readlines()
    for i in range(0,len(param_types)):
     for line in lines:
       if line.startswith("%s  "%(str(param_types[i].zfill(3)))):
          #print(line)
          for j in range(0,4):
             dih_ref[i,j] = line.split()[j+1]
             
          break 
    out=open("OPLS_ref","w+")
    out.write("OPLS refernce data\n")
    for i in range(0,len(param_types)):
      out.write("%3i %3.3f %3.3f %3.3f %3.3f\n"%(int(param_types[i]),float(dih_ref[i,0]),float(dih_ref[i,1]),float(dih_ref[i,2]),float(dih_ref[i,3])))
    #print(str(df_type[0]).zfill(3))
    #print(df_type)
    #print(dih_ref)
    return dih_ref
###################################################################################################
def dih_prep(molecule_name): 
    G1, V1, A1, D1, Q1 =dih_tag(molecule_name)
    zmat=open("%s.z"%(molecule_name),"r")
    lines=zmat.readlines()
    line=lines[dihnum].split()
    C1=line[0]  
    C2=line[4]
    C3=line[6]
    C4=line[8]
    f=open("%s.z"%(molecule_name),"r")
    lines=f.readlines()
    f.close()
    dihedrals=open("Enviroment_dihedral_list.dat","r")
    d_lines=dihedrals.readlines()
    out=open("dihzmat","w+")
    i=0
    for (n, line) in enumerate(lines):
          if V1 < n < A1:
            i=i+1 
            for l in d_lines[i-1:i]:
             di1=l.split()[0]
             di2=l.split()[1]
             di3=l.split()[2]
             di4=l.split()[3]
             d_tag=l.split()[4]
             if C1 == di1 and C2 == di2 and C3 == di3 and C4 == di4:
              Main_d=d_tag
              continue
             else:
                out.write(line[0:4]+"% 4i% 4i"%(int(d_tag),int(d_tag))+line[12:])  
          elif D1 > n > A1:
             i=i+1
             for l in d_lines[i-1:i]:
              di1=l.split()[0]
              di2=l.split()[1]
              di3=l.split()[2]
              di4=l.split()[3]
              d_tag=l.split()[4]
              out.write("% 4i% 4i% 4i% 4i% 4i% 4i\n"%(int(di1),int(di2), int(di3), int(di4), int(d_tag), int(d_tag)))       
          elif n == int(dihnum):
               out.write(line[:56]+'%11.6f'%(0)+line[67:])
          elif n == D1:
             out.write("% 4s% 4s% 4s% 4s% 4i% 4i\n                    Domain Definitions follow     (4I4)\n"%(C1,C2,C3,C4,int(Main_d),int(Main_d)))
          elif n > Q1 and BA==1 and 'X ' not in line:
             out.write(line.lower()) 
          else:
             out.write(line)
    return None  
###################################################################################################  
def new_zmat(molecule_name): 
    G1, V1, A1, D1, Q1 =dih_tag(molecule_name)
    zmat=open("%s.z"%(molecule_name),"r")
    lines=zmat.readlines()
    line=lines[dihnum].split()
    C1=line[0]   
    C2=line[4]
    C3=line[6]
    C4=line[8]
    f=open("%s.z"%(molecule_name),"r")
    lines=f.readlines()
    f.close()
    dihedrals=open("Enviroment_dihedral_list.dat","r")
    d_lines=dihedrals.readlines()
    out=open("zmat.new","w+")
    i=0
    for (n, line) in enumerate(lines):
          if V1 < n < A1:
            i=i+1 
            for l in d_lines[i-1:i]:
             d_tag=l.split()[4]
             out.write(line[0:4]+"% 4i% 4i"%(int(d_tag),int(d_tag))+line[12:]) 
          elif D1 > n > A1:
             i=i+1
             for l in d_lines[i-1:i]:
              di1=l.split()[0]
              di2=l.split()[1]
              di3=l.split()[2]
              di4=l.split()[3]
              d_tag=l.split()[4]
              out.write("% 4i% 4i% 4i% 4i% 4i% 4i\n"%(int(di1),int(di2), int(di3), int(di4), int(d_tag), int(d_tag))) 
          elif n > Q1 and 'X ' not in line:
                 out.write(line.lower()) 
          else:
             out.write(line)
    return None  
###################################################################################################    
def find_torsions_types(molecule_name, dihnum):
    global tor_list, N_torsions, param_types 
    data = open('%s.z'%(molecule_name)).read()
    count_dum = data.count('DUM')
    tor_list=[]
    G1, V1, A1, D1, Q1 =dih_tag(molecule_name)
    var_dihs=np.zeros((A1-(V1+1),2))
    i=0
    with open("%s.z"%(molecule_name),"r") as f:
     lines=f.readlines()
     for  line in lines[V1+1:A1]:
            var_dihs[i,0]=line.split()[0]
            var_dihs[i,1]=line.split()[1]
            i=i+1            
    out=open("dihedral_list.dat","w+")
    i=0 
    with open("%s.z"%(molecule_name),"r") as f:
      lines=f.readlines()
      for  line in lines[int(var_dihs[0,0]):int(var_dihs[len(var_dihs)-1,0])+1]:
        out.write("%4s %4s %4s %4s %4i\n" %(line.split()[0],line.split()[4],line.split()[6],line.split()[8],int(var_dihs[i,1])))
        i=i+1  
    out.close()
#now just extract the premade list of additional dihedrals and add them to the variabile list
    d_list=open("dihedral_list.dat","a")
    with open("%s.z"%(molecule_name),"r") as f:
     lines=f.readlines()
     for  line in lines[A1+1:D1]:
         d_list.write("%4s %4s %4s %4s %4s\n" %(line.split()[0],line.split()[1],line.split()[2],line.split()[3],line.split()[4]))     
    d_list.close()
    zmat=open("%s.z"%(molecule_name),"r")
    Z_lines=zmat.readlines()
    QM_names=[] #get the opls atom type names from the zmat
    for line in Z_lines[Q1+1:]:
      if ' X ' not in line:
        try: 
          QM_names.append(int(line.split()[0]))
          if len(str(line.split()[2])) ==1:
              name=str(line.split()[2])+' '
              QM_names.append(name.lower())
          else: 
              QM_names.append(str(line.split()[2]).lower())
        except:
            continue
    QM_names=np.reshape(QM_names, (int(len(QM_names)/2),2))
    line=Z_lines[dihnum].split()
    C1=line[0]  
    C2=line[4]
    C3=line[6]
    C4=line[8]
    d_list=open("dihedral_list.dat","r")
    D_lines=d_list.readlines()
    for line in D_lines:
     A1=line.split()[0]
     A2=line.split()[1]
     A3=line.split()[2]
     A4=line.split()[3]
     d_tag=line.split()[4] 
     if int(d_tag) not in improper_list and int(d_tag) <600:
      if C2 == A2 and C3 == A3: 
        #if d_tag not in df_type:
          tor_list.append(d_tag) 
          tor_list.append(A1)
          tor_list.append(A2)
          tor_list.append(A3)
          tor_list.append(A4)
      elif C2 == A3 and C3 == A2:
           tor_list.append(d_tag) 
           tor_list.append(A4)
           tor_list.append(A3)
           tor_list.append(A2)
           tor_list.append(A1)
    tor_list=np.reshape(tor_list, (int(len(tor_list)/5),5))
    N_torsions=len(tor_list)
    print('''| Total amount of dihedrals describing the torsion %s                                             |''' %(len(tor_list)))
    print('''| The amount of parameter types %s                                                                |''' %(len(set(tor_list[0:,0]))))
    print('''| Torsion list:                                                                                  |''')
    print('''| Type    i    j    k    l                                                                       |''')  
    for i in range(0,len(tor_list)):
     print('''| %4i %4i %4i %4i %4i                                                                       |''' %(int(tor_list[i,0]), int(tor_list[i,1]), int(tor_list[i,2]), int(tor_list[i,3]),int(tor_list[i,4])))  
    print('''| Dihedrals with the same type and different enviroment will be split and re-numbered            |''')
    #print(N_torsions)
    #found the torsions in the dihedral now test the enviroments?
    #print(count_dum)
    Atoms=G1-1-count_dum
    #print(Atoms)
    QM_tag=[]
    i=0
    with open("%s.z"%(molecule_name),"r") as f:
     lines=f.readlines()
     for  line in lines[int(Q1+2):int(Q1+2+Atoms)]:
         QM_tag.extend([line.split()[2]])         
         i=i+1
    #print(QM_tag)
    enviroment_listA=[]
    enviroment_listB=[]
    param_types=np.unique(tor_list[0:,0])
    #print(param_types)  #the amount of unique param types
    split=[]     #array of param types that appear more than once that may need spliting 
    for i in range(0,len(param_types)):
       if np.sum(tor_list[0:,0] == param_types[i]) != 1:
        split.append(param_types[i])
    #print(split)
    for i in range(0,len(split)):
        for j in range(0,len(tor_list)):
         if tor_list[j,0] == split[i]:
            if tor_list[j,1] not in enviroment_listA:
               enviroment_listA.append(tor_list[j,1])
            if tor_list[j,4] not in enviroment_listB:
               enviroment_listB.append(tor_list[j,4])
    #print(enviroment_listA)   
    #print(enviroment_listB)         
    enviroment_tagsA=[]
    enviroment_tagsB=[]
    for i in range(0,len(enviroment_listA)):
        near_list=[]
        for line in D_lines:
         if int(line.split()[4]) not in improper_list:
          if enviroment_listA[i] == line.split()[0] and line.split()[1] not in near_list and int(line.split()[1]) > count_dum:
                near_list.append(line.split()[1])
          elif enviroment_listA[i] == line.split()[1] and line.split()[2] not in near_list and int(line.split()[2]) > count_dum:
                near_list.append(line.split()[2])
          elif enviroment_listA[i] == line.split()[2] and line.split()[3] not in near_list and int(line.split()[3]) > count_dum:
                near_list.append(line.split()[3])
          elif enviroment_listA[i] == line.split()[3] and line.split()[2] not in near_list and int(line.split()[2]) > count_dum:
                near_list.append(line.split()[2])      
        if len(near_list) != 4:
         for i in range(0,4-len(near_list)):
          near_list.append(0)            
        #print(near_list)
        for i in range(0,4):  #can only have a max of three atoms bonded 
         if int(near_list[i]) !=0:
           enviroment_tagsA.append(QM_tag[int(near_list[i])-1-count_dum])
         else:
           enviroment_tagsA.append('DU')
    #print(enviroment_tagsA)
   
    
    for i in range(0,len(enviroment_listB)):
        near_list=[]
        for line in D_lines:
         if int(line.split()[4]) not in improper_list:
          if enviroment_listB[i] == line.split()[0] and line.split()[1] not in near_list and int(line.split()[1]) > count_dum:
                near_list.append(line.split()[1])
          elif enviroment_listB[i] == line.split()[1] and line.split()[2] not in near_list and int(line.split()[2]) > count_dum:
                near_list.append(line.split()[2])
          elif enviroment_listB[i] == line.split()[2] and line.split()[3] not in near_list and int(line.split()[3]) > count_dum:
                near_list.append(line.split()[3])
          elif enviroment_listB[i] == line.split()[3] and line.split()[2] not in near_list and int(line.split()[2]) > count_dum:
                near_list.append(line.split()[2])       
        if len(near_list) != 4:
         for i in range(0,4-len(near_list)):
          near_list.append(0)            
        #print(near_list)
        for i in range(0,4):  #can only have a max of four atoms bonded 
         if int(near_list[i]) !=0:
           enviroment_tagsB.append(QM_tag[int(near_list[i])-1-count_dum])
         else:
           enviroment_tagsB.append('DU')
    #print(enviroment_tagsB)
    #work out how many terms are needed for each split member of the list
    if len(split) == 0:
      print('''| All dihedral types in the torsion are different no spliting needed!                            |''')
    for i in range(0,len(split)):
      for j in range(0,len(tor_list)):
        if split[i] == tor_list[j,0]:
             posA=enviroment_listA.index(tor_list[j,1])*4
             #print(posA)
             posB=enviroment_listB.index(tor_list[j,4])*4
             #print(posB)
             envir_T1_A=enviroment_tagsA[posA:posA+4]
             envir_T1_B=enviroment_tagsB[posB:posB+4]
             #print(envir_T1_A)
             #print(envir_T1_B)
             #print(j)
             for l in range(j,len(tor_list)):
               if l != j and split[i] == tor_list[l,0]:
                 posA2=enviroment_listA.index(tor_list[l,1])*4 
                 posB2=enviroment_listB.index(tor_list[l,4])*4
                 envir_T2_A=enviroment_tagsA[posA2:posA2+4]
                 #print(envir_T2_A)
                 envir_T2_B=enviroment_tagsB[posB2:posB2+4]
                 #print(envir_T2_B)
                 compare_A=collections.Counter(envir_T1_A) == collections.Counter(envir_T2_A)
                 compare_B=collections.Counter(envir_T1_B) == collections.Counter(envir_T2_B)
                 if str(compare_A) == 'True':
                   if str(compare_B) == 'True':
                       continue
                      #print('%s and %s are the same type no split needed'%(tor_list[j], tor_list[l]))
                   else:
                       #print('%s and %s have the same type number but are not the same and need spliting'%(tor_list[j], tor_list[l]) ) 
                       #print('New tor list')
                        tor_list[l,0] = 'S%s'%(j)
                       #print(tor_list[l])
                        if 'S%s'%(j) not in split:
                            split.append('S%s'%(j))
                 else:             
                    #print('%s and %s have the same type number but are not the same and need spliting'%(tor_list[j], tor_list[l]) ) 
                    #print('New tor list')
                       tor_list[l,0] = 'S%s'%(j)
                    #print(tor_list[l])
                       if 'S%s'%(j) not in split:
                           split.append('S%s'%(j))
    #print(tor_list)
    new_tors=[] 
    for i in range(0,len(param_types)):
      if param_types[i] not in new_tors:
          new_tors.append(param_types[i])        
    for i in range(0,len(split)):
        if split[i] not in new_tors:
          new_tors.append(split[i])
          #print(new_tors)
          #print(tor_list)
    #now rewite the dihedral file with new param numbers
    out=open('Enviroment_dihedral_list.dat','w+')
    for line in D_lines:
      A1=line.split()[0]
      A2=line.split()[1]
      A3=line.split()[2]
      A4=line.split()[3]
      j=0
      for i in range(0,len(tor_list)):
       S_tag=tor_list[i,0]
       C1=tor_list[i,1]
       C2=tor_list[i,2]
       C3=tor_list[i,3]
       C4=tor_list[i,4]
       pos_S=new_tors.index(S_tag)
       dihnum=new_dihnum+pos_S  
       if A1 == C1 and A2 == C2 and A3 == C3 and A4 == C4:
         out.write("%4s %4s %4s %4s %4i\n" %(A1,A2,A3,A4,dihnum))
       elif A1==C4 and A2==C3 and A3==C2 and A4==C1:
         out.write("%4s %4s %4s %4s %4i\n" %(A1,A2,A3,A4,dihnum))
       else:
         j=j+1
         if j == len(tor_list):
          out.write(line)
      
    N_torsions=len(new_tors)    #not the same as the amount of torsions in the diehdral the amount of terms that will be refit 
    print('''| Amount of dihedral types found to be fit %s                                                     |'''%(N_torsions))
    tor_string=[]
    #print(tor_list)
    #print(param_types)
    for i in range(len(param_types)):
        for j in range(len(tor_list)):
            if str(param_types[i]) == str(tor_list[j,0]):
               tor_string.extend([QM_names[int(tor_list[j,1])-3,1], QM_names[int(tor_list[j,2])-3,1], QM_names[int(tor_list[j,3])-3,1], QM_names[int(tor_list[j,4])-3,1]])
               break
    splits=[]
    for j in range(len(tor_list)):
        if 'S' in str(tor_list[j,0]) and str(tor_list[j,0]) not in splits:
                 tor_string.extend(['XS','PL', 'IT','!X'])
                 splits.append(str(tor_list[j,0]))
    tor_string=np.reshape(tor_string, (int(len(tor_string)/4),4))
    #print(tor_string)
    return tor_list, N_torsions, param_types, tor_string
###################################################################################################
def BONDS_fit():
    print('Bonds will be fit using matlab script make sure MATLAB and Gaussian are loaded') 
    for filename in os.listdir("."):
      if filename == 'zmat.log':
        print('zmat.log found checking...')
        log=open('zmat.log','r').read()
        error=log.count('Error ')
        if error > 0:
          sys.exit("Error found in zmat.log file please check the optimisation has converged")
        else:
          print('No errors found now fitting new bonds and angles using MATLAB')
          if "lig.fchk" not in os.listdir("."):
              os.system('formchk lig.chk lig.fchk')
          os.system('cp $QuBeKit/matlab/* .')
          os.system('cp $BOSSdir/oplsaa.sb .')
          os.system("sed -i 's/..\/bonds/..\/BONDS/g' test.m")
          os.system("sed -i 's/..\/bonds/..\/BONDS/g' test.m")
          os.system("sed -i 's/vibrational_scaling= 0.957/vibrational_scaling= %s/g' test.m"%(vib_scaling))
          os.system('matlab -nojvm -nodisplay -nosplash < test.m > output.txt > /dev/null 2>&1')
          os.system('sleep 10')
          os.system('touch Modified_Scaled_Seminario.sb')
          f=open("Modified_Scaled_Seminario.sb","r")
          for (i, line) in enumerate(f):
            if "line above must be blank" in line:
              L1=i-1
          f.close()
          f=open("oplsaa.sb","r")
          for (i, line) in enumerate(f):
              if "line above must be blank" in line:
                  L2=i
          f.close()
          f=open("Modified_Scaled_Seminario.sb","r")
          lines=f.readlines()
          x=(len(lines))-2
          out=open("QuBe.sb","w+")
          f2=open("oplsaa.sb","r")
          lines2=f2.readlines()
          for (n, line) in enumerate(lines[0:x]):
                    if n < L1:
                      out.write(line[0:8].lower()+line[8:])
                    elif n == L1+1:
                      for line in lines2[1:L2]:
                        out.write(line)
                      out.write("********                        line above must be blank\r\n")
                    elif  L1+1< n < x:
                        out.write(line[0:8].lower()+line[8:])
          for line in lines2[L2+1:]:
                          out.write(line)
          print('New bonds wrote to QuBe.sb') 
          print('Now extracting xyz onetep input')
          with open('zmat.log') as log:
             for (i,line) in enumerate(log):
                if ' Standard orientation:' in line:
                   ST=i
          log=open('zmat.log')
          lines=log.readlines()
          out=open('temp.xyz','w+')
          for (i,line) in enumerate(lines):
              if i >= ST+5:
                if '-------' in line:
                  break
                else:
                  if int(line.split()[1]) == 1:
                       tag='H'
                  elif int(line.split()[1]) == 6:
                       tag='C'
                  elif int(line.split()[1]) == 7:
                       tag='N'
                  elif int(line.split()[1]) == 8:
                       tag='O'  
                  elif int(line.split()[1]) == 9:
                       tag='F'
                  elif int(line.split()[1]) == 15:
                       tag='P'
                  elif int(line.split()[1]) == 16:
                       tag='S'
                  elif int(line.split()[1]) == 17:
                       tag='Cl'
                  elif int(line.split()[1]) == 35:
                       tag='Br'
                  elif int(line.split()[1]) == 5:
                       tag='B'
                  out.write('%4s %12.6f%12.6f%12.6f\n'%(tag, float(line.split()[3]), float(line.split()[4]), float(line.split()[5])))
          out.close()        
          os.system("wc -l temp.xyz | awk '{print $1}' >> %s.xyz"%(molecule_name))
          os.system("echo '' >> %s.xyz"%(molecule_name))
          os.system('cat temp.xyz >> %s.xyz'%(molecule_name))
          os.system('rm temp.xyz')                           
    G1, V1, A1, D1, Q1 = dih_tag(molecule_name)
    zmat=open('%s.z'%(molecule_name),'r')
    out=open('%s_BA.z'%(molecule_name),'w+')
    lines=zmat.readlines()
    for (n,line) in enumerate(lines):
        if n > Q1+1:
           out.write(line.lower())
        else:
           out.write(line)
    print('New zmat (%s_BA.z) made with new bonds and angles, to use new parameters make sure QuBe.sb is present'%(molecule_name))       
###################################################################################################
def dih_tag(molecule_name):
    global G1, V1, A1, D1, Q1
    with open('%s.z'%(molecule_name)) as zmat:
     for (i, line) in enumerate(zmat):
      if "Geometry Variations follow" in line:
        G1 = i 
      elif "Variable Dihedrals follow" in line:
        V1 = i  
      elif "Additional Dihedrals follow" in line:
        A1 = i
      elif "Domain Definitions follow" in line:
        D1 = i  
      elif "QM" in line:
        Q1 = i
     return G1, V1, A1, D1, Q1 
###################################################################################################
def labels():
    global FF, NET, D1, B1, A1, DM1, NB, END
    #os.system('$BOSSdir/scripts/xSPM %s'%(zmat_name))
    with open('out','r') as sp:
     for (i, line) in enumerate(sp):
      if "OPLS Force Field Parameters" in line:
        FF = i 
      elif "Net Charge" in line:
        NET = i  
      elif "Dihedral                             Fourier Coefficients" in line:
        D1 = i  
      elif "Bond Stretching Parameters" in line:
        B1 = i
      elif "Angle Bending Parameters" in line:
        A1 = i  
      elif "Dipole Moment" in line:
        DM1 = i
      elif "Non-bonded Pairs List" in line:
         NB = i 
      elif "Checking for potential problems in the solute Z-matrix file" in line:
         END = i 
     return FF, NET, D1, B1, A1, DM1, NB, END
###################################################################################################
def dih_tag(zmat_name):  #tags from the zmat
    global G1z, V1z, A1z, D1z, Q1z
    with open('%s'%(zmat_name)+'.z') as zmat:
     for (i, line) in enumerate(zmat):
      if "Geometry Variations follow" in line:
        G1z = i 
      elif "Variable Dihedrals follow" in line:
        V1z = i  
      elif "Additional Dihedrals follow" in line:
        A1z = i
      elif "Domain Definitions follow" in line:
        D1z = i  
      elif "QM" in line:
        Q1z = i
     return G1z, V1z, A1z, D1z, Q1z
###################################################################################################
def GMX_itp():  #make the GMX itp file to go with the gro file 
    #print('Making ITP file')
    global QM_pos, QM_params
    G1z, V1z, A1z, D1z, Q1z =dih_tag(zmat_name)
    var_dihs=np.zeros((A1z-(V1z+1),2))
    i=0
    with open("%s"%(zmat_name)+'.z',"r") as f:
     lines=f.readlines()
     for  line in lines[V1z+1:A1z]:
            var_dihs[i,0]=line.split()[0]
            var_dihs[i,1]=line.split()[1]
            i=i+1
    out=open("dihedral_list.dat","w+")
    i=0 
    with open("%s"%(zmat_name)+'.z',"r") as f:
      lines=f.readlines()
      for  line in lines[int(var_dihs[0,0]):int(var_dihs[len(var_dihs)-1,0])+1]:
        out.write("%4s %4s %4s %4s %4i\n" %(line.split()[0],line.split()[4],line.split()[6],line.split()[8],int(var_dihs[i,1])))
        i=i+1  
    out.close()
#now just extract the premade list of additional dihedrals and add them to the variabile list
    d_list=open("dihedral_list.dat","a")
    with open("%s"%(zmat_name)+'.z',"r") as f:
     lines=f.readlines()
     for  line in lines[A1z+1:D1z]:
         d_list.write("%4s %4s %4s %4s %4s\n" %(line.split()[0],line.split()[1],line.split()[2],line.split()[3],line.split()[4]))     
    d_list.close()
    out=open('%s'%(resname)+'.itp','w+')
    out.write(' \n')
    out.write(';\n')
    out.write('; NEW %s itp FILE\n'%(resname))
    out.write(';\n')
    out.write('[ atomtypes ]\n')
    SP_in=open('out','r')
    lines=SP_in.readlines()
    #make array of zmat position and qm number  and elemnt type
    QM=[]
    for line in lines[0:FF]:
      for j in range(0,Atoms):
        if Atom_names[j,0] in line:
         QM.append(line.split()[0])
         QM.append(line.split()[2])
    QM=np.reshape(QM, (Atoms,2))
    QM_pos=np.concatenate((Atom_names, QM), axis=1)  
    elements=[]
    remove_digits = str.maketrans('', '', digits)
    for i in range(0,Atoms):
     element=Atom_names[i,0].translate(remove_digits)
     if len(element) == 2:
       if str(element) == 'Cl':
          elements.append(element)
       elif str(element) == 'Br':
           elements.append(element)
       else:
          elements.append(element[0])
     else:
          elements.append(element)     
             
    elements=np.reshape(elements, (Atoms,1))
    QM_pos=np.concatenate((QM_pos, elements), axis=1)
    #print(QM_pos)
    #find the qm values in the same order as Atom_names format charge | sigma | epsilon  in gro/openmm units
    QM_params=[]
    for line in lines[FF+1:NET-1]:
     for j in range(0,Atoms):
      if QM_pos[j,1] == line.split()[1]:
        QM_params.append(line.split()[2])
        QM_params.append(float(line.split()[3])/10) #sigma A to nm 
        QM_params.append(float(line.split()[4])*4.184) #epsilon kcal to kj
    QM_params=np.reshape(QM_params, (Atoms,3))
    #print(QM_params) 
    Masses=[]   
    for i in range(0,Atoms):
      if QM_pos[i,3] == 'S' :
        mass=32.0600
        Masses.append(mass)
      elif QM_pos[i,3] == 'C':
        mass=12.0110
        Masses.append(mass)
      elif QM_pos[i,3] == 'H':
        mass=1.0080
        Masses.append(mass)
      elif QM_pos[i,3] == 'F':
        mass=18.9984
        Masses.append(mass)
      elif QM_pos[i,3] == 'Cl':
        mass=35.4500
        Masses.append(mass)
      elif QM_pos[i,3] == 'N':
        mass=14.0070
        Masses.append(mass)
      elif QM_pos[i,3] == 'O':
        mass=15.9990
        Masses.append(mass)
      elif QM_pos[i,3] == 'Br':
        mass=79.9040
        Masses.append(mass)
      elif QM_pos[i,3] == 'P':
        mass=30.9738
        Masses.append(mass)
      elif QM_pos[i,3] == 'I':
        mass=126.9045
        Masses.append(mass)  
      elif QM_pos[i,3] == 'B':
        mass=10.810  
        Masses.append(mass)
      out.write('  opls_%s  %s%s%11.4f%10.3f    A    %11.5E   %10.5E\n'%(QM_pos[i,2], QM_pos[i,3], QM_pos[i,2], mass,0, float(QM_params[i,1]), float(QM_params[i,2]) ))
    #print(Masses)
    out.write('[ moleculetype ]\n')
    out.write('; Name              nrexcl\n')
    out.write('%s                  3\n'%(resname))
    out.write('[ atoms ]\n')
    out.write(';   nr       type  resnr residue  atom   cgnr     charge       mass  \n')  
    for i in range(0,Atoms):
      out.write('%6i   opls_%s      1    %s   %s      1%11.4f%11.4f\n'%(int(QM_pos[i,1])-2, QM_pos[i,2], resname, QM_pos[i,0], float(QM_params[i,0]), Masses[i]))  
    out.write('[ bonds ]\n')
    for line in lines[B1+2:A1-1]:
     out.write('%5i%6i%6i%12.4f %10.3f\n'%(int(line.split()[0])-2, int(line.split()[1])-2,1 , float(line.split()[2])/10, float(line.split()[3])*836.8))
    out.write('\n')
    out.write('[ angles ]\n')
    out.write(';  ai    aj    ak funct            c0            c1            c2            c3\n')     
    for line in lines[A1+2:DM1-1]:
      if float(line.split()[4]) != 0:
       out.write('%5i%6i%6i%6i%11.3f%11.3f\n'%(int(line.split()[0])-2, int(line.split()[1])-2, int(line.split()[2])-2, 1, float(line.split()[3]), float(line.split()[4])*8.368 ))   
    out.write('\n')
    out.write('[ dihedrals ]\n')   
    out.write('; IMPROPER DIHEDRAL ANGLES\n')
    out.write(';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n')
    i=0
    dihedrals=open('dihedral_list.dat','r')
    dih_lines=dihedrals.readlines()
    for line in lines[D1+2:B1-1]:
      if int(dih_lines[i].split()[4]) in XML_GMX_list:
       if int(dih_lines[i].split()[4]) != 0:
         F1=float(line.split()[4])*4.184
         F2=float(line.split()[5])*4.184
         F3=float(line.split()[6])*4.184
         F4=float(line.split()[7])*4.184
         C0=F2+(F1+F3)*0.5
         C1=0.5*(-F1+3*F3)
         C2=-F2+4*F4
         C3=-2*F3
         C4=-4*F4
         C5=0
         out.write('%5i%5i%5i%5i%9i%12.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n'%(int(dih_lines[i].split()[0])-2, int(dih_lines[i].split()[1])-2, int(dih_lines[i].split()[2])-2, int(dih_lines[i].split()[3])-2, 4, C0, C1, C2, C3, C4, C5))
         i=i+1
      else:
         i=i+1
    out.write('\n')
    out.write('[ dihedrals ]\n')
    out.write('; PROPER DIHEDRAL ANGLES\n')
    out.write(';  ai    aj    ak    al funct            c0            c1            c2            c3            c4            c5\n') 
    i=0   
    for line in lines[D1+2:B1-1]:
      if int(dih_lines[i].split()[4]) not in XML_GMX_list:
         F1=float(line.split()[4])*4.184
         F2=float(line.split()[5])*4.184
         F3=float(line.split()[6])*4.184
         F4=float(line.split()[7])*4.184
         C0=F2+(F1+F3)*0.5
         C1=0.5*(-F1+3*F3)
         C2=-F2+4*F4
         C3=-2*F3
         C4=-4*F4
         C5=0
         out.write('%5i%5i%5i%5i%9i%12.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n'%(int(dih_lines[i].split()[0])-2, int(dih_lines[i].split()[1])-2, int(dih_lines[i].split()[2])-2, int(dih_lines[i].split()[3])-2, 3, C0, C1, C2, C3, C4, C5))
         i=i+1
      else:
         i=i+1
    out.write('\n')
    out.write('[ pairs ]\n')
    i=0
    A=0
    for line in lines[D1+2:B1-1]:
      if int(dih_lines[i].split()[4]) not in XML_GMX_list: 
        if i!=0:
          for x in range(0,i):
              if int(dih_lines[i].split()[0])-2 == int(dih_lines[x].split()[0])-2 and int(dih_lines[i].split()[3])-2 == int(dih_lines[x].split()[3])-2:
                  A=1
              elif int(dih_lines[i].split()[0])-2 == int(dih_lines[x].split()[3])-2 and int(dih_lines[i].split()[3])-2 == int(dih_lines[x].split()[0])-2:
                   A=1
          if A==0:  
                   out.write('%5i%5i%5i\n'%(int(dih_lines[i].split()[0])-2, int(dih_lines[i].split()[3])-2, 1))          
                   i=i+1
          else: 
               i=i+1
               A=0       
        else:
             out.write('%5i%5i%5i\n'%(int(dih_lines[i].split()[0])-2, int(dih_lines[i].split()[3])-2, 1))          
             i=i+1        
      else:
           i=i+1
    out.write('\n')    
    return QM_pos, QM_params
###################################################################################################
def ITP_addsites():
    #first find all the virtual site 1-4 interaction prom the pdb file 
    Excep_pairs=[]
    constraint_no=[]
    constraint_dist=[]
    sites=0
    graph={}
    out=open('zmat.z','w+')
    zmat=open('%s.z'%(molecule_name),'r')
    lines=zmat.readlines()
    for line in lines:
      if ' X' not in line:
         out.write(line)
    out.close() 
    os.system('$BOSSdir/scripts/xZPDB zmat > /dev/null 2>&1')
    pdb=open('zmat.pdb','r')
    lines=pdb.readlines()
    for line in lines:
       if 'CONECT' in line: #get conections without extra site
         graph[int(line.split()[1])] = [int(line.split()[2])]
         for i in range(3,len(line.split())):
           if int(line.split()[i]) != 0 :
             graph[int(line.split()[1])].append(int(line.split()[i]))
    pdb.close()
    #os.system('rm zmat*')
    zmat=open('%s.z'%(molecule_name),'r')
    zmat_lines=zmat.readlines()
    for line in zmat_lines:
        if ' Geometry Variations follow' in line:
             break
        elif ' X' in line:
            site_no=int(line.split()[0])-2
            parent_no=int(line.split()[4])-2
            constraint_no.append(site_no)
            constraint_no.append(parent_no)
            constraint_dist.append(float(line.split()[5]))
            graph[site_no] = [parent_no]
            graph[parent_no].append(site_no)
          
            #print(graph) 
            for i in range(1,site_no):
                k=find_shortest_path(graph, site_no, i, path=[])
        #print(k)
       #print(len(k))
                if len(k) == 4:
                    Excep_pairs.append(i)
                    Excep_pairs.append(site_no)
    #print(Excep_pairs)
    constraint_no=np.reshape(constraint_no , (int(len(constraint_no)/2),2))
    #print(constraint_no)
    #print(constraint_dist)
    Excep_pairs=np.reshape(Excep_pairs , (int(len(Excep_pairs)/2),2))
    itp=open('%s.itp'%(resname),'r')
    lines=itp.readlines()
    out=open('%s_sites.itp'%(resname),'w+')
    for line in lines:
       if '[ moleculetype ]' in line:
          for i in range(len(QM_v_sites)):
             if i < 9:
               num=i+1
               num=str(num).zfill(2)
             else:
                num=i+1
             out.write('  opls_%s  X%s %11.4f%10.3f    A    %11.5E   %10.5E\n' %(QM_v_sites[i], num, 0.0000,0,1,0)) 
          out.write('[ moleculetype ]\n')
       elif '[ bonds ]' in line:
             for i in range(len(QM_v_sites)):
                 if i < 9:
                    num=i+1
                    num=str(num).zfill(2)
                 else:
                     num=i+1
                 out.write('%6i   opls_%s      1    %s   X%s      1%11.4f%11.4f\n'%(int(str(QM_v_sites[i])[1:])+1, QM_v_sites[i], resname,num, v_site_charges[i], 0)) 
             out.write('[ bonds ] \n')
       elif '[ pairs ]' in line: 
            out.write('[ constraints ] \n')
            out.write(';  ai  aj  funct  length_A \n')
            for i in range(len(constraint_no)):
                out.write(' %4i%4i   1     %4.3f \n'%(constraint_no[i,0], constraint_no[i,1], constraint_dist[i]/10))
            out.write('\n')
            out.write('[ pairs ]\n')
            for i in range(len(Excep_pairs)):
                out.write('%5i%5i%5i\n'%(Excep_pairs[i,0], Excep_pairs[i,1],1)) 
       else:
         out.write(line)
    out.close()  
    #print(QM_v_sites)
    #print(v_site_charges)
    
###################################################################################################
def find_shortest_path(graph, start, end, path=[]): #set up a graph for atoms in a molecule to find the virtual site pairs for itp files 
        path = path + [start]
        if start == end:
            return path
        if start not in graph.keys():
            return None
        shortest = None
        for node in graph[start]:
            if node not in path:
                newpath = find_shortest_path(graph, node, end, path)
                if newpath:
                    if not shortest or len(newpath) < len(shortest):
                        shortest = newpath
        return shortest 
###################################################################################################
def GMX_gro(pdb_name):  #make the GMX gro file from the PDB file 
    #print('Making GRO file')
    #first check to see if a virtual site version is needed
    Atoms_sites=0
    Atom_names_sites=0
    if 'extra.pdb' in os.listdir('.'):
        print("Extra sites will be added to gro and itp files but the site type must be added to the itp file before use!")
        out=open('%s_sites'%(resname)+'.gro','w+')
        out.write('NEW %s GRO FILE\n'%(resname))
        pdb=open('extra.pdb','r').read()
        Atoms_sites=pdb.count('HETATM')
        #print(Atoms)
        out.write('%5i\n'%(Atoms_sites))
        pdb=open('extra.pdb','r')
        lines=pdb.readlines()
        Atom_names_sites=[]
        for line in lines:
           if 'HETATM' in line:
               out.write('    1%s    %s   %2i%8.3f%8.3f%8.3f\n'%(resname, line.split()[2], int(line.split()[1]), float(line.split()[6])/10, float(line.split()[7])/10, float(line.split()[8])/10))
               Atom_names_sites.append(line.split()[2]) 
               o_tag=line.split()[3]
        out.write('%10.5f%10.5f%10.5f\n'%(1,1,1))
        out.close()

    out=open('%s'%(resname)+'.gro','w+')
    out.write('NEW %s GRO FILE\n'%(resname))
    pdb=open('%s'%(pdb_name),'r').read()
    Atoms=pdb.count('HETATM')
    #print(Atoms)
    out.write('%5i\n'%(Atoms))
    pdb=open('%s'%(pdb_name),'r')
    lines=pdb.readlines()
    Atom_names=[]
    for line in lines:
      if 'HETATM' in line:
       out.write('    1%s    %s   %2i%8.3f%8.3f%8.3f\n'%(resname, line.split()[2], int(line.split()[1]), float(line.split()[6])/10, float(line.split()[7])/10, float(line.split()[8])/10))
       Atom_names.append(line.split()[2]) 
       o_tag=line.split()[3]
    out.write('%10.5f%10.5f%10.5f\n'%(1,1,1))
    Atom_names=np.reshape(Atom_names, (len(Atom_names),1))
    #print(Atom_names)      
    return Atoms, Atom_names, o_tag, Atoms_sites, Atom_names_sites
###################################################################################################
def XML():   #make the xml file for openMM  can use the itp format as a lot is the same
    #print('Making XML file')
    SP_in=open('out','r')
    lines=SP_in.readlines()
    out=open('%s'%(resname)+'.xml','w+')
    out.write('<ForceField>\n')
    out.write('<AtomTypes>\n')
    #print(QM_pos)
    #print(QM_params)
    for i in range(0,Atoms):
      if QM_pos[i,3] == 'S' :
        mass=32.0600
      elif QM_pos[i,3] == 'C':
        mass=12.0110
      elif QM_pos[i,3] == 'H':
        mass=1.0080
      elif QM_pos[i,3] == 'F':
        mass=18.9984
      elif QM_pos[i,3] == 'Cl':
        mass=35.4500
      elif QM_pos[i,3] == 'N':
        mass=14.0070
      elif QM_pos[i,3] == 'O':
        mass=15.9990
      elif QM_pos[i,3] == 'Br':
        mass=79.9040
        QM_pos[i,3]='B'
      elif QM_pos[i,3] == 'P':
        mass=30.9738
      elif QM_pos[i,3] == 'I':
        mass=126.9045  
      elif QM_pos[i,3] == 'B':
        mass=10.810
      out.write('<Type name="opls_%s" class="%s%s" element="%s" mass="%8.6f" />\n'%(QM_pos[i,2], QM_pos[i,3], QM_pos[i,2], QM_pos[i,3], float(mass)))
    out.write('</AtomTypes>\n')
    out.write('<Residues>\n')
    out.write('<Residue name="%s">\n'%(resname))
    for i in range(0,Atoms):
      out.write('<Atom name="%s" type="opls_%s" />\n'%(QM_pos[i,0], QM_pos[i,2])) 
    for line in lines[B1+2:A1-1]:
      out.write('<Bond from="%s" to="%s"/>\n'%(int(line.split()[0])-3, int(line.split()[1])-3)) 
    out.write('</Residue>\n')
    out.write('</Residues>\n')
    out.write('<HarmonicBondForce>\n')
    for line in lines[B1+2:A1-1]:
      atom1=int(line.split()[0])-3 #convert to new atom number -2 for dummys and -1 for python index
      atom2=int(line.split()[1])-3 
      out.write('<Bond class1="%s%s" class2="%s%s" length="%8.6f" k="%13.6f"/>\n'%(QM_pos[atom1,3], QM_pos[atom1,2], QM_pos[atom2,3], QM_pos[atom2,2], float(line.split()[2])/10, float(line.split()[3])*836.8))
    out.write('</HarmonicBondForce>\n')
    out.write('<HarmonicAngleForce>\n')
    for line in lines[A1+2:DM1-1]:
      if float(line.split()[4]) != 0:
       atom1=int(line.split()[0])-3
       atom2=int(line.split()[1])-3
       atom3=int(line.split()[2])-3   
       out.write('<Angle class1="%s%s" class2="%s%s" class3="%s%s" angle="%8.6f" k="%10.6f"/>\n'%(QM_pos[atom1,3], QM_pos[atom1,2], QM_pos[atom2,3], QM_pos[atom2,2], QM_pos[atom3,3], QM_pos[atom3,2], float(line.split()[3])*(math.pi/180), float(line.split()[4])*8.368)) 
    out.write('</HarmonicAngleForce>\n')
    out.write('<PeriodicTorsionForce>\n')
    dihedrals=open('dihedral_list.dat','r')
    dih_lines=dihedrals.readlines()
    i=0
    for line in lines[D1+2:B1-1]:
       if int(dih_lines[i].split()[4]) != 0:
         F1=float(line.split()[4])*2.092
         F2=float(line.split()[5])*2.092
         F3=float(line.split()[6])*2.092
         F4=float(line.split()[7])*2.092
         atom1=int(dih_lines[i].split()[0])-3
         atom2=int(dih_lines[i].split()[1])-3
         atom3=int(dih_lines[i].split()[2])-3
         atom4=int(dih_lines[i].split()[3])-3
         if int(dih_lines[i].split()[4]) not in XML_GMX_list:
           out.write('<Proper class1="%s%s" class2="%s%s" class3="%s%s" class4="%s%s" k1="%8.6f" k2="%8.6f" k3="%8.6f" k4="%8.6f" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0.00" phase2="3.141592653589793" phase3="0.00" phase4="3.141592653589793"/>\n'%(QM_pos[atom1,3], QM_pos[atom1,2], QM_pos[atom2,3], QM_pos[atom2,2], QM_pos[atom3,3], QM_pos[atom3,2], QM_pos[atom4,3], QM_pos[atom4,2], F1, F2, F3, F4))
           i=i+1
         else:
            out.write('<Improper class1="%s%s" class2="%s%s" class3="%s%s" class4="%s%s" k1="%8.6f" k2="%8.6f" k3="%8.6f" k4="%8.6f" periodicity1="1" periodicity2="2" periodicity3="3" periodicity4="4" phase1="0.00" phase2="3.141592653589793" phase3="0.00" phase4="3.141592653589793"/>\n'%(QM_pos[atom1,3], QM_pos[atom1,2], QM_pos[atom2,3], QM_pos[atom2,2], QM_pos[atom3,3], QM_pos[atom3,2], QM_pos[atom4,3], QM_pos[atom4,2], F1, F2, F3, F4))
            i=i+1
       else:
         i=i+1 
    out.write('</PeriodicTorsionForce>\n')
    out.write('<NonbondedForce coulomb14scale="0.5" lj14scale="0.5">\n')
    for j in range(0,Atoms):
      out.write('<Atom type="opls_%s" charge="%8.6f" sigma="%8.6f" epsilon="%8.6f" />\n'%(QM_pos[j,2], float(QM_params[j,0]), float(QM_params[j,1]), float(QM_params[j,2]))) 
    out.write('</NonbondedForce>\n')
    out.write('</ForceField>\n')
    out.write('\n')
    out.write('\n')  
###################################################################################################
def BONDS_sub():
    out=open('bonds_script.sub','w+')
    out.write('''#!/bin/bash

#SBATCH -A DCCADD
#SBATCH --ntasks=%s

module load gaussian/09

g09 < zmat.com > zmat.log


module load MATLAB/2017a

module load Anaconda3/5.0.1 

QuBeKit.py -f bonds -t fit -z %s.z 
    '''%(processors, molecule_name))
###################################################################################################
def DIHEDRALS_sub(): 
    out=open('dihedrals_script.sub','w+')
    out.write('''#!/bin/bash
#SBATCH -A DCCADD
#SBATCH --ntasks=%s

module load gaussian/09

g09 < zmat.com > zmat.log 

    '''%(processors))
###################################################################################################
def Param_search(molecule_name, new_dihnum):  
    global BA      
    print('Looking for QuBe bonds files')
    if "QuBe.sb" in os.listdir("../../QuBe_PARAMS"):
        print('Bonds and Angles file found will be used')
        os.system('cp ../../QuBe_PARAMS/QuBe.sb .')
        BA=1
    else:
        print('No new bonds and angles found using opls values') 
        os.system('cp $BOSSdir/oplsaa.sb .')
        BA=0
    if "%s_NB_BAD.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS") or "%s_NBV_BAD.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
        use=input('Other dihedral fitted parameters found (this is the case when fitting multiple torsions) would you like to use these?:\n> ')
        if use=='yes' or use == 'y':
           if "%s_NBV_BAD.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
               print("Non-bonded parameters with virtual sites found")
               molecule_name=molecule_name+'_NBV_BAD'
           else:
               print("Non-bonded parameters found")
               molecule_name=molecule_name+'_NB_BAD'
           os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
           os.system('cp ../../QuBe_PARAMS/QuBe.par oplsaa.par')
           os.system('cp oplsaa.par oplsstore')
           par=open('oplsstore','r')
           lines=par.readlines()
           for (i,line) in enumerate(lines):
               if 'NEWQuBe' in line:
                   new_dihnum=int(line.split()[0])+1            
        else:
            os.system('cp $BOSSdir/oplsaa.par .')
            os.system('cp oplsaa.par oplsstore') 
            new_dihnum=new_dihnum
            if "%s_NBV_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                print("Non-bonded parameters with virtual sites found and will be used")
                molecule_name=molecule_name+'_NBV_BA'
                os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
            elif "%s_NB_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                print("Non-bonded parameters found and will be used")
                molecule_name=molecule_name+'_NB_BA'
                os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
            elif "%s_NBV.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                 print("Non-bonded parameters with virtual sites found and will be used")
                 molecule_name=molecule_name+'_NBV'
                 os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))    
            elif "%s_NB.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                 print("Non-bonded parameters found and will be used")
                 molecule_name=molecule_name+'_NB'
                 os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
    elif "%s_BAD.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          use=input('Other dihedral fitted parameters found (this is the case when fitting multiple torsions) would you like to use these?:\n> ')
          if use=='yes' or use == 'y':
             os.system('cp ../../QuBe_PARAMS/QuBe.par oplsaa.par')
             os.system('cp oplsaa.par oplsstore')
             par=open('oplsstore','r')
             lines=par.readlines()
             for (i,line) in enumerate(lines):
                 if 'NEWQuBe' in line:
                     new_dihnum=int(line.split()[0])+1              
          else:
               os.system('cp $BOSSdir/oplsaa.par .')
               os.system('cp oplsaa.par oplsstore')
               new_dihnum=new_dihnum 
               if "%s_NBV_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                   print("Non-bonded parameters found and will be used")
                   molecule_name=molecule_name+'_NBV_BA'
                   os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
               elif "%s_NB_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                   print("Non-bonded parameters found and will be used")
                   molecule_name=molecule_name+'_NB_BA'
                   os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
               elif "%s_NBV.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                     print("Non-bonded parameters found and will be used")
                     molecule_name=molecule_name+'_NBV'
                     os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))    
               elif "%s_NB.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                     print("Non-bonded parameters found and will be used")
                     molecule_name=molecule_name+'_NB'
                     os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
    elif "%s_D.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          use=input('Other dihedral fitted parameters found (this is the case when fitting multiple torsions) would you like to use these?:\n> ')
          if use=='yes' or use == 'y':
             os.system('cp ../../QuBe_PARAMS/QuBe.par oplsaa.par')
             os.system('cp oplsaa.par oplsstore')
             par=open('oplsstore','r')
             lines=par.readlines()
             for (i,line) in enumerate(lines):
                 if 'NEWQuBe' in line:
                     new_dihnum=int(line.split()[0])+1              
          else:
               os.system('cp $BOSSdir/oplsaa.par .')
               os.system('cp oplsaa.par oplsstore')
               new_dihnum=new_dihnum 
               if "%s_NBV_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                   print("Non-bonded parameters found and will be used")
                   molecule_name=molecule_name+'_NBV_BA'
                   os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
               elif "%s_NB_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                   print("Non-bonded parameters found and will be used")
                   molecule_name=molecule_name+'_NB_BA'
                   os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))
               elif "%s_NBV.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                     print("Non-bonded parameters found and will be used")
                     molecule_name=molecule_name+'_NBV'
                     os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))    
               elif "%s_NB.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
                     print("Non-bonded parameters found and will be used")
                     molecule_name=molecule_name+'_NB'
                     os.system("cp ../../QuBe_PARAMS/%s.z ."%(molecule_name))                     
    elif "%s_NBV_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          print('Non-bonded parameters with virtual sites found and will be used')
          molecule_name=molecule_name+'_NBV_BA'
          os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
          if "QuBe.par" in os.listdir("../../QuBe_PARAMS"): 
              print('QuBe replaced zmat and par files found')
              os.system('cp ../../QuBe_PARAMS/QuBe.par oplsaa.par')
          else:
               os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore')                        
    elif "%s_NB_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          print('Non-bonded parameters found and will be used')
          molecule_name=molecule_name+'_NB_BA'
          os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
          if "QuBe.par" in os.listdir("../../QuBe_PARAMS"): 
              print('QuBe replaced zmat and par files found')
              os.system('cp ../../QuBe_PARAMS/QuBe.par oplsaa.par')
          else:
               os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore') 
    elif "%s_NBV.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          print('only non bonded parameters found and will be used')
          molecule_name=molecule_name+'_NBV'
          os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
          os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore')       
    elif "%s_NB.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          print('only non bonded parameters found and will be used')
          molecule_name=molecule_name+'_NB'
          os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
          os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore') 
    elif "%s_BA.z"%(molecule_name) in os.listdir("../../QuBe_PARAMS"):
          molecule_name=molecule_name+'_BA'
          os.system('cp ../../QuBe_PARAMS/%s.z .'%(molecule_name))
          os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore')
    else:
          print('No other parameters found using opls')
          os.system('cp ../%s.z .'%(molecule_name))
          os.system('cp $BOSSdir/oplsaa.par .')
          os.system('cp oplsaa.par oplsstore')     
          os.system("sed -i 's/$BOSSdir\/oplsaa.par/.\/oplsaa.par/g' dihcmd")
    if BA==1:
             os.system("sed -i 's/$boxes\/oplsaa.sb/QuBe.sb/g' dihcmd")
    os.system('''sed -i 's/set lambda = "30.000 0.000 0.000"/set lambda = "%6.3f 0.000 0.000"/g' dihcmd'''%(increment))
    os.system(''' sed -i 's/$BOSSdir\/oplsaa.par/.\/oplsaa.par/g' dihcmd''')
    if "ligandqm" in os.listdir(".") or "ligandmm" in os.listdir("."):
        os.system('rm ligand*')
    if "results.dat" not in os.listdir("."):
        sys.exit("QM results missing")
    return BA, molecule_name, new_dihnum
###################################################################################################
def xyztozmat():
        Z_mat=molecule_name+'_NB.z'
        angle_pos=[0, 0, 0]
        cross_pos=[0, 0, 0]
        v_site=0
        V_sites=[]
        V_index=[]
        cross_no=9999
        out=open('no_extras.xyz','w+')
        xyz=open('xyz_with_extra_point_charges.xyz','r')
        xyz_lines=xyz.readlines()
        #find the corodinates of the heavy atom of the v-sites
        for line in xyz_lines[2:]:
            if 'X' not in line:
                out.write(line)
        out.close()
        no_v=open('no_extras.xyz','r')
        no_v_lines=no_v.readlines()
        for (n,line) in enumerate(xyz_lines[2:]):
           if 'X' in line:
                v_pos=np.array((float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ))
                v_site=v_site+1
                site=n
                V_index.append(site)
                for i in range(site):
                   if 'X' != str(xyz_lines[site+1-i].split()[0]):
                     #print(xyz_lines[site+1-i])
                     heavy_pos=np.array(( float(xyz_lines[site+1-i].split()[1]), float(xyz_lines[site+1-i].split()[2]), float(xyz_lines[site+1-i].split()[3]) ))
                     for (x,line)  in enumerate(no_v_lines):
                         if heavy_pos[0]==float(line.split()[1]) and heavy_pos[1] == float(line.split()[2]) and heavy_pos[2] == float(line.split()[3]):
                            heavy_no=x
                            #print(heavy_no)
                     bond_dist = np.linalg.norm(heavy_pos-v_pos)
                     #print(bond_dist)
                     angle_dist=2
                     for (i, line) in enumerate(no_v_lines):
                       if 'H' not in line:                                                                                             
                         angle_pos=np.array(( float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ))
                         #print(angle_pos)
                         if angle_pos[0] != heavy_pos[0] or angle_pos[1] != heavy_pos[1] or angle_pos[2] != heavy_pos[2]:
                            dist = np.linalg.norm(heavy_pos - angle_pos)
                            #print(dist)
                            if dist < angle_dist:
                               angle_dist = dist
                               angle_no=i
                               #print(angle_no)
                     angle_pos=np.array(( float(no_v_lines[angle_no].split()[1]), float(no_v_lines[angle_no].split()[2]), float(no_v_lines[angle_no].split()[3]) ))
                     ba=v_pos-heavy_pos
                     bc=angle_pos-heavy_pos
                     cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
                     angle = np.arccos(cosine_angle)
                     #print("v site parent angle atom angle %s"%(angle*180/math.pi))
                     cross_dist=3
                     for (i,line) in enumerate(no_v_lines):
                                     if i != angle_no and heavy_no != i:
                                        cross_pos=np.array(( float(line.split()[1]), float(line.split()[2]), float(line.split()[3]) ))
                                        dist = np.linalg.norm(heavy_pos - cross_pos)
                                        if dist < cross_dist:
                                           cross_dist=dist
                                           cross_no=i
                     if cross_no != 9999:
                            cross_pos=  np.array(( float(no_v_lines[cross_no].split()[1]), float(no_v_lines[cross_no].split()[2]), float(no_v_lines[cross_no].split()[3]) ))
                     #print("cross atom position %s"%cross_pos)
                     p0=v_pos
                     p1=heavy_pos
                     p2=angle_pos
                     p3=cross_pos
                     b0= -1.0*(p1-p0)
                     b1 = p2-p1
                     b2=p3-p2
                     b1 /= np.linalg.norm(b1)
                     v = b0 - np.dot(b0, b1)*b1
                     w = b2 - np.dot(b2, b1)*b1
                     x = np.dot(v, w)
                     y = np.dot(np.cross(b1, v), w)
                     torsion2=np.degrees(np.arctan2(y, x))
                     #print("torsion angle %s"%torsion2)
                     #print("torsion atom index not zmat index %s"%cross_no)
                     V_sites.append(heavy_no+3)
                     V_sites.append(bond_dist)
                     V_sites.append(angle_no+3)
                     V_sites.append(angle*180/math.pi)
                     V_sites.append(cross_no+3)
                     if abs(torsion2) > 0.001:
                           V_sites.append(torsion2)
                     else:
                           V_sites.append(0)
                     break
        V_sites=np.reshape(V_sites, (v_site,6))
        print(V_sites)
        with open(Z_mat) as zmat:
             for (i, line) in enumerate(zmat):
              if "Geometry Variations follow" in line:
                 G1z = i 
              elif "Variable Dihedrals follow" in line:
                 V1z = i  
              elif "Additional Dihedrals follow" in line:
                 A1z = i
              elif "Domain Definitions follow" in line:
                 D1z = i  
              elif "QM" in line:
                 Q1z = i
	        
        atom=0        
        zmat=open(Z_mat,'r')
        zmat_lines=zmat.readlines()
        out=open('%s_NBV.z'%(molecule_name),'w+')
        for (n,line) in enumerate(zmat_lines):
            if n==G1z-1:
               out.write(line)
               opls=int(line.split()[2])
               atom=int(line.split()[0])
            elif n==G1z:
                 for j in range(v_site):
                     out.write("%4i X0%1i %4i %4i %4i %11.6f%4i %11.6f%4i %11.6f UNK    1\n"%(atom+j+1, j+1, opls+j+1, opls+j+1, V_sites[j,0], V_sites[j,1], V_sites[j,2],  V_sites[j,3], V_sites[j,4], V_sites[j,5])) 
                 out.write(line)
            elif Q1z+1 < n < Q1z+atom:
                 try:
                    QM=int(line.split()[0][1:])
                 except:
                     pass
                 #print(QM)
                 for i in range(len(V_sites)): #make sure charges are symetric 
                     if QM == int(V_sites[i,0])-3: #test if parent if so print new charge minus the sites
                        out.write(line[:11]+"%10.6f "%(float(no_v_lines[QM].split()[4]))+line[22:])
                        break
                     elif i == len(V_sites)-1:
                         out.write(line)
            elif n == Q1z+atom:
                 for j in range(v_site):
                     out.write("%4i%3i X  %10.6f%10.6f%10.6f\n"%(opls+j+1, 0, float(xyz_lines[V_index[j]+2].split()[4]), 0, 0))         
            else:
                 out.write(line)
        out.write('\n')	
###################################################################################################
def xml_addextras():
    #weights should be left as they are 
     w1o,w2o,w3o =  1.0,0.0,0.0 # SUM SHOULD BE 1
     w1x,w2x,w3x = -1.0,1.0,0.0 # SUM SHOULD BE 0
     w1y,w2y,w3y = -1.0,0.0,1.0 # SUM SHOULD BE 0 
     Z_name=molecule_name+'.z'
     PDB_name='extra.pdb'
     v_sites=0  
     Index=[]
     QM_v_sites=[]
     v_site_charges=[]
     zmat=open(Z_name,'r')
     lines=zmat.readlines()
     for line in lines:
         if 'X0' in line:
            v_sites=v_sites+1
            Index.append(int(line.split()[4]))
            Index.append(int(line.split()[6]))
            Index.append(int(line.split()[8]))
            QM_v_sites.append(int(line.split()[2]))
     Index=np.reshape(Index, (v_sites,3))
     O_name=[]
     B_name=[]
     C_name=[]
     for x in range(len(Index)):
         Index_A=Index[x,0]
         Index_B=Index[x,1]
         Index_C=Index[x,2]
         Index_sites=np.zeros((v_sites,1))
         name_sites= []
         i=0
         bot=False
         for (n, line) in enumerate(lines):
             if n == Index_A:
                O_name.append(line.split()[1])
             elif n == Index_B:
    	          B_name.append(line.split()[1])
             elif n == Index_C:
    	          C_name.append(line.split()[1])
             elif "X0" in line:
                  name_sites.extend([line.split()[1]])
                  Index_sites[i,0]=n
                  i=i+1
             elif "QM" in line:
                  bot=True     
             elif bot:
                 try:
                   if int(line.split()[0]) in QM_v_sites:
                      v_site_charges.append(float(line.split()[3]))
                 except:
                      pass     
     #print(O_name)
     #print(B_name)
     #print(C_name)
     pos_v=np.zeros((v_sites,3))
     pos_a=np.zeros((v_sites,3))
     pos_b=np.zeros((v_sites,3))
     pos_c=np.zeros((v_sites,3)) 
     pdb=open("%s"%(PDB_name),"r")
     lines=pdb.readlines()
     for x in range(v_sites):
       for line in lines:
         if O_name[x] in line:
            for i in range(0,3):
                pos_a[x,i]=line.split()[6+i]
         elif B_name[x] in line:
              for i in range(0,3):
                  pos_b[x,i]=line.split()[6+i] 
         elif C_name[x] in line:
              for i in range(0,3):
                  pos_c[x,i]=line.split()[6+i]
         elif "X0" in line:
               for i in range(0,v_sites):
                    if line.split()[2]==name_sites[i]:
                       for j in range(0,3):
                           pos_v[i,j]=line.split()[6+j] 
     #print(pos_v)
     #print(pos_a)
     #print(pos_b)
     #print(pos_c)                      
     orig=[]
     AB=[]
     AC=[]
     Zdir=[]
     Xdir=[]
     Ydir=[]                      
     for x in range(v_sites):
        orig.append(w1o*pos_a[x]+w2o*pos_b[x]+pos_c[x]*w3o) 
        AB.append(w1x*pos_a[x]+w2x*pos_b[x]+w3x*pos_c[x]) #rb-ra
        AC.append(w1y*pos_a[x]+w2y*pos_b[x]+w3y*pos_c[x]) #rb-ra
        Zdir.append(np.cross(AB[x],AC[x]))
        #print(Zdir)
        Zdir[x] = Zdir[x]/np.sqrt(np.dot(Zdir[x],Zdir[x].reshape(3,1)))
        #print(Zdir)
        Xdir.append(AB[x]/np.sqrt(np.dot(AB[x],AB[x].reshape(3,1))))
        #print(Xdir)
        Ydir.append(np.cross(Zdir[x],Xdir[x])) 
        #print(Ydir) 
     #print(Zdir)
     #print(Ydir)
     #print(Xdir)        
     print("%s Extra sites found"%(v_sites))
     xml=open('%s.xml'%(resname),'r')
     lines=xml.readlines()
     bond=False
     out=open('%s_extra.xml'%(resname),'w+')
     for line in lines:
          if "</AtomTypes>" in line:
             for j in range(v_sites):
                 out.write('''<Type name="v-site%s"  class="X%s"   mass="0" />\n'''%(j+1,j+1))  
             out.write(line)         
          elif "<Bond from=" in line:
              if not bond:
                 for j in range(v_sites):
                   out.write('''<Atom name="X%s"  type="v-site%s" />\n'''%(j+1,j+1))          
                 bond=True
                 out.write(line) 
              else:
                  out.write(line)         
          elif "</NonbondedForce>" in line:
             for j in range(v_sites):
                 out.write('''<Atom type="v-site%s"  charge="%9.6f" sigma="1.000000" epsilon="0.000000" />\n'''%(j+1, v_site_charges[j]))
             out.write(line)
          elif "</Residue>" in line:                        
                for i in range(v_sites):
                    p1 = np.dot((pos_v[i]-orig[i]),Xdir[i].reshape(3,1))
                    p2 = np.dot((pos_v[i]-orig[i]),Ydir[i].reshape(3,1))
                    p3 = np.dot((pos_v[i]-orig[i]),Zdir[i].reshape(3,1))
                    out.write( '''<VirtualSite type="localCoords" index="%s" atom1="%s" atom2="%s" atom3="%s" wo1="%s" wo2="%s" wo3="%s" wx1="%s" wx2="%s" wx3="%s" wy1="%s" wy2="%s" wy3="%s" p1="%4.3f" p2="%4.3f" p3="%4.3f"/>\n'''%(int(Index_sites[i]-3),Index[i,0]-3,Index[i,1]-3,Index[i,2]-3,w1o,w2o,w3o,w1x,w2x,w3x,w1y,w2y,w3y,p1/10,p2/10,p3/10))
                out.write(line)
          else:
              out.write(line)  
     out.close()  
     return QM_v_sites, v_site_charges
###################################################################################################
def FREQ(): #same method as in the Modified Seminario paper
    where=os.getcwd()
    if 'XML+GMX' in where or 'QuBe_PARAMS' in where:
        print("Performing virbational frequency QM MM comparison")
        print("Looking for zmat.log file")
        if 'zmat.log' in os.listdir("."):
           print("Found")
        elif 'XML+GMX' in where:
             if 'zmat.log' in os.listdir("../../BONDS/"):
                 print("Found")
                 os.system("cp ../../BONDS/zmat.log .")
        elif 'zmat.log' in os.listdit("../BONDS/"):
                 print("Found")
                 os.system("cp ../BONDS/zmat.log .")
        #find the frequncies in the QM file less the 6 vibrational modes
        log=open('zmat.log','r')
        lines=log.readlines()
        QMfreq=[]
        for line in lines:
            if 'Frequencies' in line: 
               for i in range(3):
                 try:
                   if float(line.split()[i+2]) >= 0:
                      QMfreq.append(float(line.split()[i+2]))
                 except:
                          pass
        print("QM Frequencies")
        print(QMfreq)  
        #Now find the MM frequencies using parameters acording to the file name
        os.system("cp $BOSSdir/scripts/FREQcmd .")
        if  "_BAD" in molecule_name:  #assuming you always have new bonds as charges have little effect alone
            os.system("sed -i 's/$BOSSdir\/oplsaa.par/QuBe.par/g' FREQcmd")
        if "_BA" in molecule_name:
            os.system("sed -i 's/$BOSSdir\/oplsaa.sb/QuBe.sb/g' FREQcmd")
        os.system("cp %s.z optzmat"%(molecule_name))
        os.system("csh FREQcmd > /dev/null 2>&1")
        MMin=open('out','r')
        lines=MMin.readlines()
        MMfreq=[]
        dash=0
        Mode=False
        for line in lines: 
          if '   Mode ' in line:
              Mode=True
              dash=0
          elif '-----' in line:
              dash=dash+1
              if dash==2:
                 Mode=False
          elif Mode:
             MMfreq.append(float(line.split()[0]))
        print("MM Frequencies")
        print(MMfreq)
        print(len(QMfreq))
        print(len(MMfreq))
        #Now calculate the errors in the frequencies
        assert len(QMfreq)==len(MMfreq), "The QM and MM calculations must produce the same amount of frquencies"
        N=(len(QMfreq)+6)/3 #find the amount of atoms from the QM frequencies 
        front=1/(3*N-6)
        QMfreq=np.array(QMfreq)
        bot=vib_scaling*QMfreq
        top=bot-MMfreq
        mean_percent=100*front*sum(abs(top/bot))
        mean_error=front*sum(abs(top))
        print("Mean percentage error across frequencies = %5.3f"%(mean_percent))
        print("Mean unsigned error across frequencies = %6.3f"%(mean_error))
    else:
        sys.exit('Please move to QuBe_PARAMS or XML+GMX folder before running')
###################################################################################################
def find_shortest_path(graph, start, end, path=[]):
        path = path + [start]
        if start == end:
            return path
        if start not in graph.keys():
            return None
        shortest = None
        for node in graph[start]:
            if node not in path:
                newpath = find_shortest_path(graph, node, end, path)
                if newpath:
                    if not shortest or len(newpath) < len(shortest):
                        shortest = newpath
        return shortest
###################################################################################################
def OPLS_LJ(system, Excep_pairs, normal_pairs):
            forces = {system.getForce(index).__class__.__name__: system.getForce(index) for index in range(system.getNumForces())}
            nonbonded_force = forces['NonbondedForce']
            lorentz = CustomNonbondedForce(
                            'epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)*4.0')
            lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
            lorentz.addPerParticleParameter('sigma')
            lorentz.addPerParticleParameter('epsilon')
            lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
            system.addForce(lorentz)
            LJset = {}
            for index in range(nonbonded_force.getNumParticles()):
                charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
                #print(nonbonded_force.getParticleParameters(index))
                LJset[index] = (sigma, epsilon)
                lorentz.addParticle([sigma, epsilon])
                nonbonded_force.setParticleParameters(
                index, charge, sigma, epsilon * 0)
            for i in range(nonbonded_force.getNumExceptions()):
                (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
                #ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
                # FORCE
                lorentz.addExclusion(p1, p2)
                if eps._value != 0.0:
                   sig14 = sqrt(LJset[p1][0] * LJset[p2][0])
                   eps14 = sqrt(LJset[p1][1] * LJset[p2][1])
                   nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
                for x in range(len(Excep_pairs)): #scale 14 interactions
                    if p1 == Excep_pairs[x,0] and p2 == Excep_pairs[x,1] or p2 == Excep_pairs[x,0] and p1 == Excep_pairs[x,1]:
                       charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                       charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                       q = charge1*charge2*0.5
                       #print('charge %s'%q)
                       sig14 = sqrt(sigma1 * sigma2)*0.5
                       eps = sqrt(epsilon1 * epsilon2)*0.5          
                       nonbonded_force.setExceptionParameters(i, p1, p2, q , sig14, eps)
                for x in range(len(normal_pairs)):
                    if p1 == normal_pairs[x,0] and p2 == normal_pairs[x,1] or p2 == normal_pairs[x,0] and p1 == normal_pairs[x,1]:
                       charge1, sigma1, epsilon1 = nonbonded_force.getParticleParameters(p1)
                       charge2, sigma2, epsilon2 = nonbonded_force.getParticleParameters(p2)
                       q = charge1*charge2
                       #print(q)
                       sig14 = sqrt(sigma1 * sigma2)
                       eps = sqrt(epsilon1 * epsilon2)          
                       nonbonded_force.setExceptionParameters(i, p1, p2, q , sig14, eps)         

            return system
###################################################################################################
def SING_P(): #single point boss openMM comparison 
    where=os.getcwd()
    if 'XML+GMX' in where:
        print('Performing BOSS and OpenMM single point energy comparison')
        if "_BAD" in molecule_name:
          os.system('cp $BOSSdir/scripts/SPMcmd .')
          os.system("sed -i 's/$BOSSdir\/oplsaa.par/QuBe.par/g' SPMcmd")
          os.system("sed -i 's/$BOSSdir\/oplsaa.sb/QuBe.sb/g' SPMcmd")
        elif "_BA" in molecule_name:
            os.system('cp $BOSSdir/scripts/SPMcmd .')
            os.system("sed -i 's/$BOSSdir\/oplsaa.sb/QuBe.sb/g' SPMcmd")
        else: 
            os.system('cp $BOSSdir/oplsaa* .')             
            os.system('cp $BOSSdir/scripts/SPMcmd .')
        os.system('csh SPMcmd > /dev/null 2>&1') #do boss single point calc
        #find the energy in kcal/mol
        boss=open('sum','r')
        lines=boss.readlines()
        boss_energy=float(lines[0].split()[5])
        print('BOSS single point energy = %10.4f kcal/mol'%boss_energy)
        #set up the openMM single point
    temperature = 298.15 * kelvin
    pdb = PDBFile('MOL.pdb')
    modeller = Modeller(pdb.topology, pdb.positions)
    if 'MOL_extra.xml' in os.listdir("."): #are virtual sites used?
        forcefield = ForceField('MOL_extra.xml')
        modeller.addExtraParticles(forcefield)
    else:
        forcefield = ForceField('MOL.xml')
    PDBFile.writeFile(modeller.topology, modeller.positions, open('MOL_modeller.pdb', 'w'))
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff,  constraints=None)
    graph={}
    site_no=0
    Excep_pairs=[]
    normal_pairs=[]
    pdb=open('MOL_modeller.pdb','r')
    lines=pdb.readlines()
    for line in lines:
        if 'CONECT' in line:
           graph[int(line.split()[1])-1] = [int(line.split()[2])-1]
           for i in range(3,len(line.split())):
               graph[int(line.split()[1])-1].append(int(line.split()[i])-1)
    pdb.close()
    zmat=open('%s.z'%molecule_name,'r')
    lines=zmat.readlines()
    for line in lines:
        if ' X0' in line:
            site_no=int(line.split()[0])-3
            graph[int(line.split()[0])-3] = [int(line.split()[4])-3] #add site and parent 
            graph[int(line.split()[4])-3].append(int(line.split()[0])-3) #add site to parent list
            #print(graph)
            for i in range(site_no):
                k=find_shortest_path(graph, site_no, i, path=[])
                if len(k) == 4:
                   Excep_pairs.append(i)
                   Excep_pairs.append(site_no)
                elif len(k) >= 5:
                   normal_pairs.append(i)
                   normal_pairs.append(site_no)
    Excep_pairs=np.reshape(Excep_pairs , (int(len(Excep_pairs)/2),2))
    #print(Excep_pairs)
    normal_pairs=np.reshape(normal_pairs , (int(len(normal_pairs)/2),2))             
    #print(normal_pairs)
    system = OPLS_LJ(system, Excep_pairs, normal_pairs)         #OPLS mixing correction
    integrator = LangevinIntegrator(
        temperature, 5 / picosecond,  0.0005 * picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    state = simulation.context.getState(getEnergy=True)
    energy = str(state.getPotentialEnergy())
    print('OpenMM single point energy = %10.4f kcal/mol'%(float(energy[:-6])/4.184))
###################################################################################################
def Replace(): #Replace QuBeKit optimised torsions
    if 'QuBe_PARAMS' not in os.getcwd():
        sys.exit('Please use replace inside of the QuBe_PARAMS folder with the highst parameterised molecule zmat this should be name_NB/V_BA.z')
    os.system('cp $QuBeKit/bin/QuBeKit_torsions.dat .')
    G1, V1, A1, D1, Q1 =dih_tag(molecule_name)
    var_dihs=np.zeros((A1-(V1+1),2))
    i=0
    f=open("%s.z"%(molecule_name),"r")    
    lines=f.readlines()
    f.close()
    for  line in lines[V1+1:A1]:  #find the amount of variable dihedrals in the molecule
            var_dihs[i,0]=line.split()[0]
            var_dihs[i,1]=line.split()[1]
            i=i+1
    dihedral_list=list() #this will store a list of the torsions using zmat numbering and the OPLS torsion type
    i=0 
    for  line in lines[int(var_dihs[0,0]):int(var_dihs[len(var_dihs)-1,0])+1]:
         dihedral_list.extend([int(line.split()[0]), int(line.split()[4]), int(line.split()[6]), int(line.split()[8]), int(var_dihs[i,1])])
         i=i+1
    #print(dihedral_list)  
#now just extract the premade list of additional dihedrals and add them to the variabile list
    for  line in lines[A1+1:D1]:
         dihedral_list.extend([int(line.split()[0]), int(line.split()[1]), int(line.split()[2]), int(line.split()[3]), int(line.split()[4]) ]) 
    dihedral_list=np.reshape(dihedral_list, (int(len(dihedral_list)/5),5))
    #print(dihedral_list)
    QM_names=[] #get the opls atom type names from the zmat
    for line in lines[Q1+1:]:
      if ' X ' not in line:
        try: 
          QM_names.append(int(line.split()[0]))
          if len(str(line.split()[2]))==1:
             name=str(line.split()[2])+' '
             QM_names.append(str(name).lower())
          else:
             QM_names.append(str(line.split()[2]).lower())
        except:
            continue
    QM_names=np.reshape(QM_names, (int(len(QM_names)/2),2))
    #print(QM_names)
    QuBe=open('QuBeKit_torsions.dat','r') #get the QuBe torsions
    QuBe_lines=QuBe.readlines()
    os.system('cp $BOSSdir/oplsaa.par QuBe.par')
    par=open('QuBe.par','a')
    found_torsions=[]
    for line in QuBe_lines[2:]:
        write_par=True
        for x in range(len(dihedral_list)):   #check to see it the molecule has torsions in the QuBe file
            if int(dihedral_list[x,4]) not in improper_list:
               if QM_names[int(dihedral_list[x,0])-3,1]+'-'+QM_names[int(dihedral_list[x,1])-3,1]+'-'+QM_names[int(dihedral_list[x,2])-3,1]+'-'+QM_names[int(dihedral_list[x,3])-3,1] in line or QM_names[int(dihedral_list[x,3])-3,1]+'-'+QM_names[int(dihedral_list[x,2])-3,1]+'-'+QM_names[int(dihedral_list[x,1])-3,1]+'-'+QM_names[int(dihedral_list[x,0])-3,1] in line:
                  #print(dihedral_list[x])
                  #print(line.split()[0])
                  dihedral_list[x,4]=int(line.split()[0]) #if it does change the type number to the QuBe value
                  if write_par: #for each torsion in the molecule append the QuBe value to a par file once!
                     par.write('%s  % 4.3f    % 4.3f    % 4.3f    % 4.3f      %s-%s-%s-%s     **NEWQuBe**\n' %(int(line.split()[0]), float(line.split()[1]), float(line.split()[2]), float(line.split()[3]), float(line.split()[4]), QM_names[int(dihedral_list[x,0])-3,1], QM_names[int(dihedral_list[x,1])-3,1], QM_names[int(dihedral_list[x,2])-3,1], QM_names[int(dihedral_list[x,3])-3,1]))
                     write_par=False
                     found_torsions.extend([QM_names[int(dihedral_list[x,0])-3,1], QM_names[int(dihedral_list[x,1])-3,1], QM_names[int(dihedral_list[x,2])-3,1], QM_names[int(dihedral_list[x,3])-3,1] ])
    #print(dihedral_list)
    par.close()
    found_torsions=np.reshape(found_torsions, (int(len(found_torsions)/4),4))
    if len(found_torsions)!=0:
       os.system('mv %s.z %s.dat'%(molecule_name, molecule_name))
       print('Torsions found and replaced:')
       for i in range(len(found_torsions)):
           print('%s-%s-%s-%s'%(found_torsions[i,0], found_torsions[i,1], found_torsions[i,2], found_torsions[i,3]))
       for x in range(len(dihedral_list)):
           if dihedral_list[x,4] not in improper_list and dihedral_list[x,4] < 600:
              
              zmat=open('%s.z'%(molecule_name),'w+')
              break
           elif x == len(dihedral_list)-1:
              zmat=open('%sD.z'%(molecule_name),'w+')
              print('All variable torsions replaced')
              
       #Now we need to write a new zmat with the replaced QuBe torsions
       i=0
       for (n, line) in enumerate(lines):
           if V1 < n < A1: 
              zmat.write(line[0:4]+"% 4i% 4i"%(int(dihedral_list[i,4]), int(dihedral_list[i,4]))+line[12:])
              i=i+1 
           elif D1 > n > A1:
               zmat.write("% 4i% 4i% 4i% 4i% 4i% 4i\n"%(int(dihedral_list[i,0]),int(dihedral_list[i,1]), int(dihedral_list[i,2]), int(dihedral_list[i,3]), int(dihedral_list[i,4]), int(dihedral_list[i,4]))) 
               i=i+1
           else:
              zmat.write(line)
       print('''New replace file has overwriten: %s.z
Old molecule data saved as %s.dat'''%(molecule_name, molecule_name))
       zmat.close()
    else:
          os.system('rm QuBe.par') 
          print('No torsions could be replaced')
    os.system('rm QuBeKit_torsions.dat')   
###################################################################################################
#command list
###################################################################################################
if args.function == 'bonds' and args.type == 'write' and args.zmat:
   print('PDB file needed to write bonds input file making PDB')
   os.system('$BOSSdir/scripts/xZPDB %s'%(molecule_name))
   where=os.getcwd()
   if 'BONDS' in where:
       BONDS_write()
       os.system('cp %s.z Zmat.z'%(molecule_name))
       if args.submission:
          BONDS_sub()
   else: 
      if 'BONDS' not in os.listdir("."):
          os.system('mkdir BONDS')
      os.system('cp %s.pdb BONDS/'%(molecule_name))
      os.system('cp %s.z BONDS/Zmat.z'%(molecule_name))
      os.chdir('BONDS')    
      BONDS_write()
      if args.submission:
         BONDS_sub()
###################################################################################################        
elif args.function == 'bonds' and args.type == 'write' and args.PDB:
     print('Zmat file is also needed for fitting stage and will be moved make sure it is present and the order matches the coresponding pdb file')
     where=os.getcwd()
     if 'BONDS' in where:
         BONDS_write()
         os.system('cp ../%s.z .'%(molecule_name))
         if args.submission:
            BONDS_sub()
     else:
        if 'BONDS' not in os.listdir("."):
            os.system('mkdir BONDS') 
        os.system('cp %s.pdb BONDS/'%(molecule_name))
        if '%s.z'%(molecule_name) not in os.listdir("."):
            sys.exit('%s.z file not found please make a matching set of zmat and pdb files'%(molecule_name))
        else:
             os.system('cp %s.z BONDS/'%(molecule_name))    
             os.chdir('BONDS')    
             BONDS_write()
             if args.submission:
                BONDS_sub()
     settings()
###################################################################################################           
elif args.function == 'bonds' and args.type == 'fit' and args.zmat:
     where=os.getcwd()
     if 'BONDS' in where:
         if "Zmat.z" not in os.listdir("."):
             if "%s.z"%(molecule_name) in os.listdir("."):
                 os.system('cp %s.z Zmat.z'%(molecule_name))
             else:
                  os.system('cp ../%s.z Zmat.z'%(molecule_name))
         BONDS_fit()
         os.system('rm Zmat.z')
         os.chdir('../')
         if 'CHARGES' not in os.listdir("."):
             os.system('mkdir CHARGES')
         os.system('mv BONDS/%s.xyz CHARGES/'%(molecule_name))
         os.system('cp %s.z CHARGES/'%(molecule_name))
         print('onetep xyz file wrote to new folder charges')
         if 'QuBe_PARAMS' not in os.listdir("."):
             os.system('mkdir QuBe_PARAMS')
         os.system('cp BONDS/QuBe.sb QuBe_PARAMS')
         os.system('mv BONDS/%s_BA.z QuBe_PARAMS/'%(molecule_name))
     else:
         os.chdir('BONDS')
         os.system('cp ../%s.z .'%(molecule_name))
         os.system('cp %s.z Zmat.z'%(molecule_name))
         BONDS_fit()
         os.system('rm Zmat.z')
         os.chdir('../')
         if 'CHARGES' not in os.listdir("."):
             os.system('mkdir CHARGES')
         os.system('mv BONDS/%s.xyz CHARGES/'%(molecule_name))
         os.system('cp %s.z CHARGES/'%(molecule_name))
         print('onetep xyz file wrote to new folder charges') 
         if 'QuBe_PARAMS' not in os.listdir("."):
             os.system('mkdir QuBe_PARAMS')
         os.system('cp BONDS/QuBe.sb QuBe_PARAMS/') 
         os.system('mv BONDS/%s_BA.z QuBe_PARAMS/'%(molecule_name))   
###################################################################################################
elif args.function == 'bonds' and args.type == 'fit' and args.PDB:
     print('Zmat is needed for fitting, it must be in the same order as that used to write the g09 input files! Run again using zmat input -z')          
###################################################################################################        
elif args.function == 'dihedrals' and args.type == 'write' and args.zmat:
   zmat_numbers=DIHEDRALS_find()
   dihnum=input('Enter the dihedral number to be scaned around, or enter all to wirte files for the minimum required torsions:\n> ')
   if dihnum=='all':
      print('Writing scan files for all variable dihedrals')
      G1, V1, A1, D1, Q1 = dih_tag(molecule_name)
      for x in range(0,len(zmat_numbers)):
          dihnum=zmat_numbers[x]
          os.system('mkdir SCAN_%s'%(zmat_numbers[x]))
          os.system('cp %s.z SCAN_%s/'%(molecule_name, zmat_numbers[x]))
          os.chdir('SCAN_%s'%(zmat_numbers[x]))
          DIHEDRALS_write()
          os.chdir('../')
   else:   
    where=os.getcwd()
    if 'SCAN_%s'%(dihnum) in where:
       G1, V1, A1, D1, Q1 = dih_tag(molecule_name)
       DIHEDRALS_write()
    else:
       if 'SCAN_%s'%(dihnum) not in os.listdir("."):
           os.system('mkdir SCAN_%s'%(dihnum))
       os.system('cp %s.z SCAN_%s/'%(molecule_name, dihnum))
       os.chdir('SCAN_%s'%(dihnum))
       G1, V1, A1, D1, Q1 = dih_tag(molecule_name)
       DIHEDRALS_write()
   settings() 
###################################################################################################
elif args.function == 'dihedrals' and args.type == 'write' and args.PDB:
     print('Zmat file is needed to write the scan input files please try again using the zmat identifier -z')
###################################################################################################
elif args.function == 'dihedrals' and args.type == 'analyse':
     where=os.getcwd()
     grid=np.linspace(dihstart, (increment*numscan)-increment, num=numscan)
     grid=np.reshape(grid, (int(len(grid)),1))
     if 'SCAN' in where:
         print('Checking zmat.log files for normal termination and collecting results')
         print(f'Progress key: {Fore.YELLOW}+{Style.RESET_ALL}=Running, {Fore.GREEN}O{Style.RESET_ALL}=Finished with no errors, {Fore.RED}X{Style.RESET_ALL}=Error found check output')
         errors, result=DIHEDRALS_analyse()
         grid=np.concatenate((grid,result), axis=1)
         for x in range(len(grid)):
             if grid[x,1]=='O':
                print(f'Angle  %5s  {Fore.GREEN}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
             elif grid[x,1]=='X':
                 print(f'Angle  %5s  {Fore.RED}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
             else:
                  print(f'Angle  %5s  {Fore.YELLOW}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
         if errors==0:
            print('Results collected and wrote to results.dat file with no errors')
         else:
            print('Errors found please check output')
     else:
          scans=[]
          for filename in sorted(os.listdir(".")):
              if "SCAN" in filename:
                  scans.append(filename)
          print('Multiple dihedral scans found collecting results for all scans:')        
          print('Checking zmat.log files for normal termination and collecting results')
          print(f'Progress key: {Fore.YELLOW}+{Style.RESET_ALL}=Running, {Fore.GREEN}O{Style.RESET_ALL}=Finished with no errors, {Fore.RED}X{Style.RESET_ALL}=Error found check output')
          for i in range(0,len(scans)):
              os.chdir('%s'%(scans[i]))
              #print(os.getcwd())
              errors, result=DIHEDRALS_analyse()
              grid=np.concatenate((grid,result), axis=1)
              os.chdir('../')
          print('            '+'%-2s ' * len(scans)%(tuple(scans)))
          for x in range(len(grid)):  
              print('Angle  %5s   '%(grid[x,0]), end=" ")
              for y in range(len(scans)):
                  if grid[x,y+1]=='X':
                     if y==len(scans)-1:
                        print(f'{Fore.RED}X{Style.RESET_ALL}')
                     else:
                         print(f'{Fore.RED}X{Style.RESET_ALL}      ', end=" ")
                  elif grid[x,y+1]=='O':
                      if y==len(scans)-1:
                          print(f'{Fore.GREEN}O{Style.RESET_ALL}')
                      else: 
                           print(f'{Fore.GREEN}O{Style.RESET_ALL}      ', end=" ")
                  else:
                     if y==len(scans)-1:
                        print(f'{Fore.YELLOW}+{Style.RESET_ALL}')
                     else:
                        print(f'{Fore.YELLOW}+{Style.RESET_ALL}      ', end=" ")      
###################################################################################################
elif args.function == 'dihedrals' and args.type == 'check' and args.zmat:
     where=os.getcwd()
     if 'SCAN' in where:
         print('Checking torsion profile')
         os.system('mkdir Checker')
         os.chdir("Checker")
         os.system('cp ../results.dat .')
         os.system('cp $BOSSdir/testjobs/dihdrive/dihcmd .')
         BA, molecule_name, new_dihnum = Param_search(molecule_name, new_dihnum)
         starting_zmat(molecule_name)
         cmd_prep(dihnum)
         QM_prep(Q_file)
         run_boss()
         torsionparams = np.zeros((4,4))
         dih_ref = np.zeros((4,4))
         sumerror, penalty =errorcal(0,dih_ref,torsionparams, T_wieght)
         plot()
         os.system("cp Plot Starting_plot")
         os.system("rm Plot ")
         clean()
         os.system("rm ligandqm dihzmat") 
###################################################################################################
elif args.function == 'dihedrals' and args.type == 'fit' and args.zmat:
     where=os.getcwd()
     grid=np.linspace(dihstart, (increment*numscan)-increment, num=numscan)
     grid=np.reshape(grid, (int(len(grid)),1))
     if "FITTING" in where:
              if "results.dat" in os.listdir("."):
                  option=input('would you like to perform fitting in this folder with the parameters that are all ready here(yes), or search for new parameters and run fitting(no)?:\n> ')
                  if option=='yes' or option=='y':
                     if "ligandqm" in os.listdir(".") or "ligandmm" in os.listdir("."):
                         os.system('rm ligand*')
                     if "QuBe.sb" in os.listdir("."):
                         BA=1
                     else:
                         BA=0
                     DIHEDRALS_fit()
                     if "_NBV_BAD" in molecule_name:
                         os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                     elif "_NB_BAD" in molecule_name:
                         os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                     elif "_BAD" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                     elif "_NBV_BA" in molecule_name:
                         os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))      
                     elif "_NB_BA" in molecule_name:
                         os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
                     elif "_BA" in molecule_name:    
                           os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
                     elif "_NBV" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))      
                     elif "_NB" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))      
                     else:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))
                     print('''==================================================================================================''')  
                  elif option=='no' or option=='n':
                       BA, molecule_name, new_dihnum = Param_search(molecule_name, new_dihnum)
                       DIHEDRALS_fit()
                       if "_NBV_BAD" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                       elif "_NB_BAD" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                       elif "_BAD" in molecule_name:
                             os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
                       elif "_NBV_BA" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))      
                       elif "_NB_BA" in molecule_name:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
                       elif "_BA" in molecule_name:    
                            os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
                       elif "_NBV" in molecule_name:
                             os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))     
                       elif "_NB" in molecule_name:
                             os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))     
                       else:
                           os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))
                       print('''==================================================================================================''')
     elif 'SCAN' in where:
         if "results.dat" not in os.listdir("."):
             print('Checking zmat.log files for normal termination and collecting results')
             print(f'Progress key: {Fore.YELLOW}+{Style.RESET_ALL}=Running, {Fore.GREEN}O{Style.RESET_ALL}=Finished with no errors, {Fore.RED}X{Style.RESET_ALL}=Error found check output')
             errors, result=DIHEDRALS_analyse()
             grid=np.concatenate((grid,result), axis=1)
             for x in range(len(grid)):
                if grid[x,1]=='O':
                   print(f'Angle  %5s  {Fore.GREEN}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
                elif grid[x,1]=='X':
                   print(f'Angle  %5s  {Fore.RED}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
                else:
                    print(f'Angle  %5s  {Fore.YELLOW}%s{Style.RESET_ALL}'%(grid[x,0], grid[x,1]))
             if errors==0:
                print('Results collected and wrote to results.dat file with no errors')
             else:
                 sys.exit('Errors found please check output')
         print('Results found')
         print('Now starting fitting')
         os.system('mkdir FITTING')
         os.system('cp results.dat %s.z FITTING/'%(molecule_name))
         os.chdir('FITTING')
         print('Using same torsion as scan for fitting %s'%(dihnum))
         os.system('cp $BOSSdir/testjobs/dihdrive/dihcmd .')
         BA, molecule_name, new_dihnum = Param_search(molecule_name, new_dihnum)
         DIHEDRALS_fit()
         if "_NBV_BAD" in molecule_name:
             os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
         elif "_NB_BAD" in molecule_name:
             os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
         elif "_BAD" in molecule_name:
               os.system('cp zmat.new ../../QuBe_PARAMS/%s.z'%(molecule_name))
         elif "_NBV_BA" in molecule_name:
             os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))      
         elif "_NB_BA" in molecule_name:
             os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
         elif "_BA" in molecule_name:    
               os.system('cp zmat.new ../../QuBe_PARAMS/%sD.z'%(molecule_name))
         elif "_NBV" in molecule_name:
               os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))      
         elif "_NB" in molecule_name:
               os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))      
         else:
              os.system('cp zmat.new ../../QuBe_PARAMS/%s_D.z'%(molecule_name))
         print('''==================================================================================================''')  
     else:
          scans=[]
          for filename in sorted(os.listdir(".")):
              if "SCAN" in filename:
                  scans.append(filename)
          print('Multiple dihedral scans found collecting results for all scans:')        
          print('Checking zmat.log files for normal termination and collecting results')
          print(f'Progress key: {Fore.YELLOW}+{Style.RESET_ALL}=Running, {Fore.GREEN}O{Style.RESET_ALL}=Finished with no errors, {Fore.RED}X{Style.RESET_ALL}=Error found check output')
          for i in range(0,len(scans)):
              os.chdir('%s'%(scans[i]))
              #print(os.getcwd())
              errors, result=DIHEDRALS_analyse()
              grid=np.concatenate((grid,result), axis=1)
              os.chdir('../')
          print('            '+'%-2s ' * len(scans)%(tuple(scans)))
          for x in range(len(grid)):  
              print('Angle  %5s   '%(grid[x,0]), end=" ")
              for y in range(len(scans)):
                  if grid[x,y+1]=='X':
                     if y==len(scans)-1:
                        print(f'{Fore.RED}X{Style.RESET_ALL}')
                     else:
                         print(f'{Fore.RED}X{Style.RESET_ALL}      ', end=" ")
                  elif grid[x,y+1]=='O':
                      if y==len(scans)-1:
                          print(f'{Fore.GREEN}O{Style.RESET_ALL}')
                      else: 
                           print(f'{Fore.GREEN}O{Style.RESET_ALL}      ', end=" ")
                  else:
                     if y==len(scans)-1:
                        print(f'{Fore.YELLOW}+{Style.RESET_ALL}')
                     else:
                        print(f'{Fore.YELLOW}+{Style.RESET_ALL}      ', end=" ")   
          print('To fit torsions run QuBeKit.py -f dihedrals -t fit -z %s.z in each scan folder in your choice of order'%(molecule_name))       
###################################################################################################
elif args.function == 'charges' and args.type == 'write' and args.zmat:
     where=os.getcwd()
     if "CHARGES" in where:
         if "%s.xyz"%(molecule_name) in os.listdir("."):
             print('%s.xyz found now writing onetep input files'%(molecule_name))
             os.system('cp /nobackup/proj/dccadd/qlj-examples/input_files/* .')
             os.system('rm AZ1.xyz')
             os.system("sed -i 's/AZ1/%s/g' run_onetep"%(molecule_name))         
         else:
              print('%s.xyz  not found make sure this has been generated from the frequency calculation'%(molecule_name))
     elif "CHARGES" in os.listdir("."):
           os.chdir('CHARGES')
           if "%s.xyz"%(molecule_name) in os.listdir("."):
             print('%s.xyz found now writing onetep input files'%(molecule_name))
             os.system('cp /nobackup/proj/dccadd/qlj-examples/input_files/* .')
             os.system('rm AZ1.xyz')
             os.system("sed -i 's/AZ1/%s/g' run_onetep"%(molecule_name))         
           else:
              print('%s.xyz  not found make sure this has been generated from the frequency calculation'%(molecule_name)) 
     else:
           print('CHARGES folder not found, make sure to run a QM frequency calculation first!')
###################################################################################################
elif args.function == 'charges' and args.type == 'fit' and args.zmat: 
     where=os.getcwd()
     if "iter_1" in where:
       if "ddec_error_message" not in os.listdir("."):
         if "%s.z"%(molecule_name) in os.listdir("."):
            os.system('cp %s.z zmat'%(molecule_name))
            os.system('cp $QuBeKit/onetep/lj_script .')
            os.system('./lj_script')
            os.system('mv zmat_ddec %s_NB.z'%(molecule_name))
            if 'xyz_with_extra_point_charges.xyz' in os.listdir("."):
                xyz=open('xyz_with_extra_point_charges.xyz').read()
                if 'X ' in xyz:
                    print('Extra sites found in xyz atempting to convert to zmat')
                    xyztozmat()
                    if 'QuBe.sb' in os.listdir("../../QuBe_PARAMS/"):
                        out_V=open('%s_NBV_BA.z'%(molecule_name),'w+')
                        G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NBV')
                        z_in=open('%s_NBV.z'%(molecule_name),'r')
                        lines=Z_in.readlines()
                        for (n ,line) in enumerate(lines):
                            if n > Q1:
                               out_V.write(line.lower())
                            else:
                                out_V.write(line)
                        out_V.close()
                        os.system('cp %s_NBV_BA.z ../../QuBe_PARAMS/'%(molecule_name))
                    os.system('cp %s_NBV.z ../../QuBe_PARAMS/'%(molecule_name))
                    print('Zmat with and without extra sites made and moved')
            out=open('%s_NB_BA.z'%(molecule_name),'w+')
            G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NB')
            z_in=open('%s_NB.z'%(molecule_name),'r')
            lines=Z_in.readlines()
            for (n ,line) in enumerate(lines):
                if n > Q1:
                    out.write(line.lower())
                else:
                    out.write(line)
            out.close()
            os.system('cp %s_NB.z ../../QuBe_PARAMS/'%(molecule_name))
            os.system('cp %s_NB_BA.z ../../QuBe_PARAMS/'%(molecule_name))
     elif "iter_1" in os.listdir("."):
           os.chdir('iter_1')
           if "ddec_error_message" not in os.listdir("."):
               if "%s.z"%(molecule_name) in os.listdir("."):
                   os.system('cp %s.z zmat'%(molecule_name))
               else:
                    os.system("cp ../../%s.z ."%(molecule_name))
                    os.system("cp %s.z zmat"%(molecule_name))
               os.system('cp $QuBeKit/onetep/lj_script .')
               os.system('./lj_script')
               os.system('mv zmat_ddec %s_NB.z'%(molecule_name))
               if 'xyz_with_extra_point_charges.xyz' in os.listdir("."):
                   xyz=open('xyz_with_extra_point_charges.xyz').read()
                   if 'X ' in xyz:
                       print('Extra sites found in xyz atempting to convert to zmat')
                       xyztozmat()
                       if 'QuBe.sb' in os.listdir("../../QuBe_PARAMS/"):
                          out_V=open('%s_NBV_BA.z'%(molecule_name),'w+')
                          G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NBV')
                          z_in=open('%s_NBV.z'%(molecule_name),'r')
                          lines=z_in.readlines()
                          for (n ,line) in enumerate(lines):
                              if n > Q1:
                                 out_V.write(line.lower())
                              else:
                                  out_V.write(line)
                          out_V.close()
                          os.system('cp %s_NBV_BA.z ../../QuBe_PARAMS/'%(molecule_name))
                       os.system('cp %s_NBV.z ../../QuBe_PARAMS/'%(molecule_name))
                       print('Zmat with and without extra sites made and moved')
               out=open('%s_NB_BA.z'%(molecule_name),'w+')
               G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NB')
               z_in=open('%s_NB.z'%(molecule_name),'r')
               lines=Z_in.readlines()
               for (n ,line) in enumerate(lines):
                   if n > Q1:
                            out.write(line.lower())
                   else:
                          out.write(line)
               out.close()
               os.system('cp %s_NB_BA.z ../../QuBe_PARAMS/'%(molecule_name))
               os.system('cp %s_NB.z ../../QuBe_PARAMS/'%(molecule_name))
     elif "CHARGES" in os.listdir("."):
           os.chdir('CHARGES/iter_1')
           if "ddec_error_message" not in os.listdir("."):
               if "%s.z"%(molecule_name) in os.listdir("."):
                   os.system('cp %s.z zmat'%(molecule_name))
               else:
                    os.system("cp ../../%s.z ."%(molecule_name))
                    os.system("cp %s.z zmat"%(molecule_name))
               os.system('cp $QuBeKit/onetep/lj_script .')
               os.system('./lj_script')
               os.system('mv zmat_ddec %s_NB.z'%(molecule_name))
               if 'xyz_with_extra_point_charges.xyz' in os.listdir("."):
                   xyz=open('xyz_with_extra_point_charges.xyz').read()
                   if 'X ' in xyz:
                       print('Extra sites found in xyz atempting to convert to zmat')
                       xyztozmat()
                       if 'QuBe.sb' in os.listdir("../../QuBe_PARAMS/"):
                           out_V=open('%s_NBV_BA.z'%(molecule_name),'w+')
                           G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NBV')
                           z_in=open('%s_NBV.z'%(molecule_name),'r')
                           lines=z_in.readlines()
                           for (n ,line) in enumerate(lines):
                               if n > Q1:
                                  out_V.write(line.lower())
                               else:
                                   out_V.write(line)
                           out_V.close()
                           os.system('cp %s_NBV_BA.z ../../QuBe_PARAMS/'%(molecule_name))
                       os.system('cp %s_NBV.z ../../QuBe_PARAMS/'%(molecule_name))
                       print('Zmat with and without extra sites made and moved')
               out=open('%s_NB_BA.z'%(molecule_name),'w+')
               G1, V1, A1, D1, Q1 = dih_tag(molecule_name+'_NB')
               z_in=open('%s_NB.z'%(molecule_name),'r')
               lines=z_in.readlines()
               for (n ,line) in enumerate(lines):
                   if n > Q1:
                            out.write(line.lower())
                   else:
                          out.write(line)
               out.close()
               os.system('cp %s_NB_BA.z ../../QuBe_PARAMS/'%(molecule_name))
               os.system('cp %s_NB.z ../../QuBe_PARAMS/'%(molecule_name))
###################################################################################################
if args.XML and args.zmat:
   print('Making openMM XML file and GMX .gro and itp files for the molecule')
   print('Make sure the zmat you want to use is entered using -z')
   print('The molecule tag name will be the first three letters of the molecule')
   where=os.getcwd()
   print(molecule_name)
   if "QuBe_PARAMS" in where:
       os.system('mkdir XML+GMX')
       if "_BAD" in molecule_name:
          os.system('cp QuBe.sb QuBe.par %s.z XML+GMX/'%(molecule_name))
          os.chdir('XML+GMX')
          os.system('cp $BOSSdir/scripts/SPMcmd .')
          os.system("sed -i 's/$BOSSdir\/oplsaa.par/QuBe.par/g' SPMcmd")
          os.system("sed -i 's/$BOSSdir\/oplsaa.sb/QuBe.sb/g' SPMcmd")
       elif "_BA" in molecule_name:
            os.system('cp QuBe.sb %s.z XML+GMX/'%(molecule_name))
            os.chdir('XML+GMX')
            os.system('cp $BOSSdir/scripts/SPMcmd .')
            os.system("sed -i 's/$BOSSdir\/oplsaa.sb/QuBe.sb/g' SPMcmd")
       else: 
            os.system('cp %s.z XML+GMX/'%(molecule_name))
            os.system('cp $BOSSdir/oplsaa* XML+GMX/')            
            os.chdir('XML+GMX')  
            os.system('cp $BOSSdir/scripts/SPMcmd .')
   else:
       sys.exit('Please use this function in the QuBe_PARAMS folder')
   os.system('$BOSSdir/scripts/xZPDB %s > /dev/null 2>&1'%(molecule_name))      
   os.system('cp %s.z optzmat'%(molecule_name))
   print('Doing Single Point calculation')
   os.system('csh SPMcmd > /dev/null 2>&1')
   zmat_name=molecule_name
   pdb_name=molecule_name+'.pdb'
   if '_NBV' in molecule_name:
      print('Extra sites found and will be wrote to xml!')
      os.system('mv %s extra.pdb'%(pdb_name))
      extra_pdb=open('extra.pdb','r')
      extra_lines=extra_pdb.readlines()
      out=open(pdb_name,'w+')
      for line in extra_lines:
           if 'X0' not in line:
               out.write(line)
      out.close()
   resname='MOL'
   FF, NET, D1, B1, A1, DM1, NB, END = labels()
   print('Making PRO file')
   Atoms, Atom_names, o_tag, Atoms_sites, Atom_names_sites =GMX_gro(pdb_name)
   print('Done')
   print('Making ITP file')
   QM_pos, QM_params = GMX_itp()
   print('Done')
   print('Making XML file')
   XML()
   print('Done')
   if '_NBV' in molecule_name:
       print('Now edditing the xml to include extra sites')
       QM_v_sites, v_site_charges=xml_addextras()
       ITP_addsites()
       out=open('Zmat.z','w+')
       zin=open('%s.z'%(molecule_name),'r')
       lines=zin.readlines()
       for line in lines:
           if 'X' not in line:
              out.write(line)
       out.close()
       os.system('$BOSSdir/scripts/xZPDB Zmat > /dev/null 2>&1')
       os.system('mv Zmat.pdb %s.pdb'%(molecule_name))
       os.system('rm Zmat.z')
   os.system("sed -i 's/%s/%s/g' %s.pdb"%(o_tag, resname, molecule_name))
   os.system('cp %s.pdb %s.pdb'%(molecule_name, resname))
###################################################################################################   
if args.frequency and args.zmat:
   FREQ()
###################################################################################################
if args.singlepoint and args.zmat:
   SING_P()   
################################################################################################### 
if args.replace and args.zmat: 
   print('Looking for matching molecule and QuBe torsions')
   Replace() 
###################################################################################################





