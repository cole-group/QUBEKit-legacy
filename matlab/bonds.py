import os 

f=open("Modified_Scaled_Seminario.sb","r")
for (i, line) in enumerate(f):
 if "line above must be blank" in line:
   L1=i-1
f.close()
#print(L1)



f=open("oplsaa.sb","r")
for (i, line) in enumerate(f):
 if "line above must be blank" in line:
   L2=i
f.close()
#print(L2)



f=open("Modified_Scaled_Seminario.sb","r")
lines=f.readlines()
x=(len(lines))-2
out=open("new.sb","w+")
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
      
