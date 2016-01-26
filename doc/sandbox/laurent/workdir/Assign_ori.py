import numpy
import matplotlib.pyplot as plt
import matplotlib
from copy import deepcopy
import random

#############################  Attribution des orientations aux sets  ##############################################################

ipath = "inp_grain.inp"  # chemin du fichier "fix-ger.inp"


inp=open(ipath,"r")
content=inp.readlines()
inp.close()

inp2=open(ipath,"r")


### Comptage du nbre de sets ###
#n=0
#while ("End Part" not in content[n]):
#      n=n+1
           
        
i=0
j=-1   
ligne=[]  
sets=[]  
l=[] 
elements=[]
size=[]
tens=file("sets.inp","w")

## Recuperation et ecriture des sets

for ligne in inp2:
    j=j+1
    k=j+1   
    if ligne[0:15]=="*ELSET, ELSET=G":    #Localisation des Elsets
       i=i+1     # Nombre de sets
       while ("*ELSET" not in content[k]):
             l=content[k]
             elements.append(l)     # creation du set par ajout des lignes delements
             k=k+1
             if ("** ELEMENT SURFACES" in content[k]):
                 break             # sortie du while quand plus de sets
       temp=','.join(map(str,elements))# conversion de la liste elements en str 
       tempo=temp.split(',')# separation du str en cellules delements
       tempp = []
       for toto in xrange(len(tempo)):
         if tempo[toto]!= '\n':
           tempp.append(tempo[toto])

       length=len(tempp)
       print length
       size.append(length)
       
       for e in range(0,length,1):
              print >> tens, '          S(%d,%d) = %s' %(i,e+1,tempp[e])  # ecriture des sets et elements
#
       elements=[]
       l=[]  
   
sizemax=max(size)     # nombre max delements dans un set
tens.seek(0,2)
print >> tens, '         RETURN'
print >> tens, '         END' 
tens.close() 
inp2.close()


# recuperation des orientations
ipath = "angle_CuBe2.txt" 
inp=open(ipath,"r")
content=inp.readlines()
random.shuffle(content)
inp.close()

ori=file("orientations.inp","w")    
l=[] 
for k in range(0,i,1):
    l=content[k]
    temp=''.join(map(str,l))
#    print temp
    tempp=temp.split( )
#    print tempp
    for e in range(0,3,1):
        print >> ori, '          Ang(%d,%d) = %s' %(k+1,e+1,tempp[e])  # ecriture des sets et elements
    l=[]

ori.seek(0,2)           #se placer en fin de fichier
print >> ori, '         RETURN'
print >> ori, '         END' 
ori.close()


# Creation de la routine de stockage des angles 
sub=open("orientations.inp","r+")
content=sub.readlines()
sub.close()
content.insert(0, '          SUBROUTINE  Angles(Ang)'+'\n'+'C'+'\n'
                  +'         IMPLICIT NONE'+ '\n'
#                  +'         INTEGER I, J, nblock, N, nori '+'\n'
                  +'         REAL*8 Ang('+ str(i) +','+ str(3) +')'+'\n'+'C'+'\n')
#                  
inp=open("Angles.for","w")
inp.writelines(content)
inp.close()    
      
      
# Creation de la routine Sets
sub=open("sets.inp","r+")
content=sub.readlines()
sub.close()

content.insert(0, '          SUBROUTINE  Sets(S)'+'\n'+'C'+'\n'
                  +'         IMPLICIT NONE'+ '\n'
#                  +'         INTEGER I, J, nblock, N, nori '+'\n'
                  +'         REAL*8 S('+ str(i) +','+ str(sizemax) +')'+'\n'+'C'+'\n')

                  
inp=open("Sets.for","w")
inp.writelines(content)
inp.close()


# Creation de la routine dattribution de orientation a lelement
sub2=open("Ori.for","w+")
content=sub2.readlines()
sub2.close()
content.insert(0, '          SUBROUTINE  Ori(nblock, N, Nel, Q)'+'\n'+'C'+'\n'
                  +'         IMPLICIT NONE'+ '\n'
                  +'         INTEGER I, J, nblock, N, Nel, nset '+'\n'
                  +'         REAL*8 S('+ str(i) +','+ str(sizemax) +'), Q(nblock,3,3), Ang('+ str(i) +','+ str(3) +'), phi1, phi, phi2'+'\n'+'C'+'\n'
                  +'         call Sets(S)           !Sets d''elements'+'\n'
                  +'         DO I=1,'+ str(i) +'\n'
                  +'              DO J=1,'+ str(sizemax) +'\n'
                  +'                    IF (Nel==S(I,J)) THEN ' +'\n'
                  +'                        nset=I      '+'\n'
                  +'                        EXIT'        +'\n'
                  +'                    ENDIF'+'\n'
                  +'              END DO'+'\n'
                  +'         END DO'+'\n'+'C'+'\n'
                  +'         call Angles(Ang)'+'\n'+'C'+'\n'
                  +'              phi1 = Ang(nset,1) '+'\n'
                  +'              phi = Ang(nset,2) '+'\n'
                  +'              phi2 = Ang(nset,3) '+'\n'+'C'+'\n'
                  +'         call Korientinit(phi1,phi,phi2,Q,nblock,N)'+'\n'+'C'+'\n'
                  +'         RETURN'+'\n'
                  +'         END'+'\n')

sub2=open("Ori.for","w")
sub2.writelines(content)
sub2.close()
