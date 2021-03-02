from pylab import *
import numpy as np
import scipy
from scipy.linalg import null_space

def Fermi(x):
        return 1.0/(exp(x)+1)

def Bose(x):
	if(x>1e-4):
        	return 1.0/(exp(x)-1)
	else:
		return 1e5

A=zeros((16,16))
c=zeros(16)

Gamma_l=0.2
Gamma_r=0.2
Gamma=0.2
Gamma_g=0.2

V_D=3.0
Vin_A=3.0
Vin_B=0.0
E_P1=Vin_A
E_P2=Vin_B
E_N1=1.5*V_D-Vin_A
E_N2=1.5*V_D-Vin_B
mu_l=0.0
mu_r=mu_l-V_D
Vout=0.0
kBT=1.0
Cg=200.0

tint=10
T=5000000
Ntot=int(T/tint)
output=np.zeros(Ntot)
current=np.zeros(Ntot)
dissipation=np.zeros(Ntot)
p_P1=np.zeros(Ntot)
p_P2=np.zeros(Ntot)
p_N1=np.zeros(Ntot)
p_N2=np.zeros(Ntot)

for i in range(Ntot):
	mu_g=0.0-Vout

	k_N2l=Gamma_l*Fermi((E_N2-mu_l)/kBT)
	k_lN2=Gamma_l*(1.0-Fermi((E_N2-mu_l)/kBT))
	k_P1r=Gamma_r*Fermi((E_P1-mu_r)/kBT)
	k_rP1=Gamma_r*(1.0-Fermi((E_P1-mu_r)/kBT))
	k_P2r=Gamma_r*Fermi((E_P2-mu_r)/kBT)
	k_rP2=Gamma_r*(1.0-Fermi((E_P2-mu_r)/kBT))
	k_P1g=Gamma_g*Fermi((E_P1-mu_g)/kBT)
	k_gP1=Gamma_g*(1.0-Fermi((E_P1-mu_g)/kBT))	
	k_P2g=Gamma_g*Fermi((E_P2-mu_g)/kBT)
	k_gP2=Gamma_g*(1.0-Fermi((E_P2-mu_g)/kBT))
	k_N1g=Gamma_g*Fermi((E_N1-mu_g)/kBT)
	k_gN1=Gamma_g*(1.0-Fermi((E_N1-mu_g)/kBT))
	if(E_N1>E_N2):
		k_N1N2=Gamma*Bose((E_N1-E_N2)/kBT)
		k_N2N1=Gamma*(1+Bose((E_N1-E_N2)/kBT))
	else:
		k_N2N1=Gamma*Bose((E_N2-E_N1)/kBT)
		k_N1N2=Gamma*(1+Bose((E_N2-E_N1)/kBT))
	if(E_N1>E_P1):
		k_N1P1=Gamma*Bose((E_N1-E_P1)/kBT)
		k_P1N1=Gamma*(1+Bose((E_N1-E_P1)/kBT))
	else:
		k_P1N1=Gamma*Bose((E_P1-E_N1)/kBT)
		k_N1P1=Gamma*(1+Bose((E_P1-E_N1)/kBT))
	if(E_N1>E_P2):
		k_N1P2=Gamma*Bose((E_N1-E_P2)/kBT)
		k_P2N1=Gamma*(1+Bose((E_N1-E_P2)/kBT))
	else:
		k_P2N1=Gamma*Bose((E_P2-E_N1)/kBT)
		k_N1P2=Gamma*(1+Bose((E_P2-E_N1)/kBT))
	if(E_P1>E_P2):
		k_P1P2=Gamma*Bose((E_P1-E_P2)/kBT)
		k_P2P1=Gamma*(1+Bose((E_P1-E_P2)/kBT))
	else:
		k_P2P1=Gamma*Bose((E_P2-E_P1)/kBT)
		k_P1P2=Gamma*(1+Bose((E_P2-E_P1)/kBT))
	
	A[1][0]=k_P1r+k_P1g
	A[2][0]=k_P2r+k_P2g
	A[3][0]=k_N1g
	A[4][0]=k_N2l
	A[0][0]=-1.0*(A[1][0]+A[2][0]+A[3][0]+A[4][0])
	A[0][1]=k_rP1+k_gP1
	A[2][1]=k_P2P1
	A[3][1]=k_N1P1
	A[5][1]=k_P2r+k_P2g
	A[6][1]=k_N1g
	A[7][1]=k_N2l
	A[1][1]=-1.0*(A[0][1]+A[2][1]+A[3][1]+A[5][1]+A[6][1]+A[7][1])
	A[0][2]=k_rP2+k_gP2
	A[1][2]=k_P1P2
	A[3][2]=k_N1P2
	A[5][2]=k_P1r+k_P1g
	A[8][2]=k_N1g
	A[9][2]=k_N2l
	A[2][2]=-1.0*(A[0][2]+A[1][2]+A[3][2]+A[5][2]+A[8][2]+A[9][2])
	A[0][3]=k_gN1
	A[1][3]=k_P1N1
	A[2][3]=k_P2N1
	A[4][3]=k_N2N1
	A[6][3]=k_P1r+k_P1g
	A[8][3]=k_P2r+k_P2g
	A[10][3]=k_N2l
	A[3][3]=-1.0*(A[0][3]+A[1][3]+A[2][3]+A[4][3]+A[6][3]+A[8][3]+A[10][3])
	A[0][4]=k_lN2
	A[3][4]=k_N1N2
	A[7][4]=k_P1r+k_P1g
	A[9][4]=k_P2r+k_P2g
	A[10][4]=k_N1g
	A[4][4]=-1.0*(A[0][4]+A[3][4]+A[7][4]+A[9][4]+A[10][4])
	A[1][5]=k_rP2+k_gP2
	A[2][5]=k_rP1+k_gP1
	A[6][5]=k_N1P2
	A[8][5]=k_N1P1
	A[13][5]=k_N2l
	A[14][5]=k_N1g
	A[5][5]=-1.0*(A[1][5]+A[2][5]+A[6][5]+A[8][5]+A[13][5]+A[14][5])
	A[1][6]=k_gN1
	A[3][6]=k_rP1+k_gP1
	A[5][6]=k_P2N1
	A[7][6]=k_N2N1
	A[8][6]=k_P2P1
	A[12][6]=k_N2l
	A[14][6]=k_P2r+k_P2g
	A[6][6]=-1.0*(A[1][6]+A[3][6]+A[5][6]+A[7][6]+A[8][6]+A[12][6]+A[14][6])
	A[1][7]=k_lN2
	A[4][7]=k_rP1+k_gP1
	A[6][7]=k_N1N2
	A[9][7]=k_P2P1
	A[10][7]=k_N1P1
	A[12][7]=k_N1g
	A[13][7]=k_P2r+k_P2g
	A[7][7]=-1.0*(A[1][7]+A[4][7]+A[6][7]+A[9][7]+A[10][7]+A[12][7]+A[13][7])
	A[2][8]=k_gN1
	# this term was missed!
	A[3][8]=k_gP2+k_rP2
	A[5][8]=k_P1N1
	A[6][8]=k_P1P2
	A[9][8]=k_N2N1
	A[11][8]=k_N2l
	A[14][8]=k_P1r+k_P1g
	A[8][8]=-1.0*(A[2][8]+A[3][8]+A[5][8]+A[6][8]+A[9][8]+A[11][8]+A[14][8])
	A[2][9]=k_lN2
	A[4][9]=k_rP2+k_gP2
	A[7][9]=k_P1P2
	A[8][9]=k_N1N2
	A[10][9]=k_N1P2
	A[11][9]=k_N1g
	A[13][9]=k_P1r+k_P1g
	A[9][9]=-1.0*(A[2][9]+A[4][9]+A[7][9]+A[8][9]+A[10][9]+A[11][9]+A[13][9])
	A[3][10]=k_lN2
	A[4][10]=k_gN1
	A[7][10]=k_P1N1
	A[9][10]=k_P2N1
	A[11][10]=k_P2r+k_P2g
	A[12][10]=k_P1r+k_P1g
	A[10][10]=-1.0*(A[3][10]+A[4][10]+A[7][10]+A[9][10]+A[11][10]+A[12][10])
	# this term was missed!
	A[8][11]=k_lN2
	A[9][11]=k_gN1
	A[10][11]=k_rP2+k_gP2
	A[12][11]=k_P1P2
	A[13][11]=k_P1N1
	A[15][11]=k_P1r+k_P1g
	A[11][11]=-1.0*(A[8][11]+A[9][11]+A[10][11]+A[12][11]+A[13][11]+A[15][11])
	# deleted a term here
	A[6][12]=k_lN2
	A[7][12]=k_gN1
	A[10][12]=k_rP1+k_gP1
	A[11][12]=k_P2P1
	A[13][12]=k_P2N1
	A[15][12]=k_P2r+k_P2g
	A[12][12]=-1.0*(A[6][12]+A[7][12]+A[10][12]+A[11][12]+A[13][12]+A[15][12])
	A[5][13]=k_lN2
	A[7][13]=k_rP2+k_gP2
	A[9][13]=k_rP1+k_gP1
	A[11][13]=k_N1P1
	A[12][13]=k_N1P2
	A[14][13]=k_N1N2
	A[15][13]=k_N1g
	A[13][13]=-1.0*(A[5][13]+A[7][13]+A[9][13]+A[11][13]+A[12][13]+A[14][13]+A[15][13])
	A[5][14]=k_gN1
	A[6][14]=k_rP2+k_gP2
	A[8][14]=k_rP1+k_gP1
	A[13][14]=k_N2N1
	A[15][14]=k_N2l
	A[14][14]=-1.0*(A[5][14]+A[6][14]+A[8][14]+A[13][14]+A[15][14])
	A[11][15]=k_rP1+k_gP1
	A[12][15]=k_rP2+k_gP2
	A[13][15]=k_gN1
	A[14][15]=k_lN2
	A[15][15]=-1.0*(A[11][15]+A[12][15]+A[13][15]+A[14][15])

	#print(np.linalg.cond(A))
	p=null_space(A)
	total=0.0
	for j in range(16):
		total+=p[j][0]
	for j in range(16):
		p[j][0]=p[j][0]*1.0/total
	p_P1[i]=p[1][0]+p[5][0]+p[6][0]+p[7][0]+p[12][0]+p[13][0]+p[14][0]+p[15][0]
	p_P2[i]=p[2][0]+p[5][0]+p[8][0]+p[9][0]+p[11][0]+p[13][0]+p[14][0]+p[15][0]
	p_N1[i]=p[3][0]+p[6][0]+p[8][0]+p[10][0]+p[11][0]+p[12][0]+p[14][0]+p[15][0]
	p_N2[i]=p[4][0]+p[7][0]+p[9][0]+p[10][0]+p[11][0]+p[12][0]+p[13][0]+p[15][0]

	J1=k_N2l*(1-p_N2[i])
	J2=k_lN2*p_N2[i]
	J3=k_P1r*(1-p_P1[i])
	J4=k_rP1*p_P1[i]
	J5=k_P2r*(1-p_P2[i])
	J6=k_rP2*p_P2[i]
	J7=k_P1g*(1-p_P1[i])
	J8=k_gP1*p_P1[i]
	J9=k_P2g*(1-p_P2[i])
	J10=k_gP2*p_P2[i]
	J11=k_N1g*(1-p_N1[i])
	J12=k_gN1*p_N1[i]

	Jg=J8-J7+J12-J11+J10-J9
	
	Vout-=1.0*Jg*tint/Cg
	output[i]=Vout
	current[i]=Jg*tint
	dissipation[i]=1.0*(J1-J2)*tint*(mu_l-mu_g)+1.0*(J3-J4+J5-J6)*tint*(mu_r-mu_g)

	print(i*tint,output[i],current[i],dissipation[i])

#file=open('result_VinA=0_VinB=0_Vd=5','w')
#for i in range(Ntot):
#        file.write("%lf %lf     %g      %g	%g	%g	%g	%g\n" %(i*tint,output[i],current[i],dissipation[i],p_P1[i],p_P2[i],p_N1[i],p_N2[i]))
#file.close()
