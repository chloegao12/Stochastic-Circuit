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

def NOT_propagation(Vin,Vout,V_D):
	A=zeros((4,4))
	c=zeros(4)

	Gamma_l=0.2
	Gamma_r=0.2
	Gamma=0.2
	Gamma_g=0.2

	E_P=Vin
	E_N=1.5*V_D-Vin
	mu_l=0.0
	mu_r=mu_l-V_D	
	kBT=1.0
	Cg=200.0

	mu_g=0.0-Vout

	k_Nl=Gamma_l*Fermi((E_N-mu_l)/kBT)
	k_lN=Gamma_l*(1.0-Fermi((E_N-mu_l)/kBT))
	k_rP=Gamma_r*(1.0-Fermi((E_P-mu_r)/kBT))
	k_Pr=Gamma_r*Fermi((E_P-mu_r)/kBT)
	k_Ng=Gamma_g*Fermi((E_N-mu_g)/kBT)
	k_gN=Gamma_g*(1.0-Fermi((E_N-mu_g)/kBT))
	k_Pg=Gamma_g*Fermi((E_P-mu_g)/kBT)
	k_gP=Gamma_g*(1.0-Fermi((E_P-mu_g)/kBT))
	if(E_N>E_P):
		k_NP=Gamma*Bose((E_N-E_P)/kBT)
		k_PN=Gamma*(1+Bose((E_N-E_P)/kBT))
	else:
		k_PN=Gamma*Bose((E_P-E_N)/kBT)
		k_NP=Gamma*(1+Bose((E_P-E_N)/kBT))

	A[1][0]=k_Pr+k_Pg
	A[2][0]=k_Nl+k_Ng
	A[0][0]=-A[1][0]-A[2][0]
	A[0][1]=k_rP+k_gP
	A[2][1]=k_NP
	A[3][1]=k_Nl+k_Ng
	A[1][1]=-A[0][1]-A[2][1]-A[3][1]
	A[0][2]=k_lN+k_gN
	A[1][2]=k_PN
	A[3][2]=k_Pr+k_Pg
	A[2][2]=-A[0][2]-A[1][2]-A[3][2]
	A[1][3]=k_lN+k_gN
	A[2][3]=k_rP+k_gP
	A[3][3]=-A[1][3]-A[2][3]

	p=null_space(A)
	sum=p[0][0]+p[1][0]+p[2][0]+p[3][0]
	p[0][0]=p[0][0]/sum
	p[1][0]=p[1][0]/sum
	p[2][0]=p[2][0]/sum
	p[3][0]=p[3][0]/sum
	p_N=p[2][0]+p[3][0]
	p_P=p[1][0]+p[3][0]
	J1=k_gN*p_N-k_Ng*(1-p_N)
	J2=k_gP*p_P-k_Pg*(1-p_P)
	J3=k_Nl*(1.0-p_N)-k_lN*p_N
	J4=k_rP*p_P-k_Pr*(1.0-p_P)

	Vout-=1.0*(J1+J2)*tint/Cg
	dissipation=1.0*J3*(mu_l-mu_g)-1.0*J4*(mu_r-mu_g)

	return Vout,dissipation

def NAND_propagation(Vin_A,Vin_B,Vout,V_D):
	A=zeros((16,16))
	c=zeros(16)

	Gamma_l=0.2
	Gamma_r=0.2
	Gamma=0.2
	Gamma_g=0.2

	mu_l=0.0
	mu_r=mu_l-V_D
	kBT=1.0
	Cg=200.0

	E_P1=Vin_A
	E_P2=Vin_B
	E_N1=1.5*V_D-Vin_A
	E_N2=1.5*V_D-Vin_B

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
	A[8][11]=k_lN2
	A[9][11]=k_gN1
	A[10][11]=k_rP2+k_gP2
	A[12][11]=k_P1P2
	A[13][11]=k_P1N1
	A[15][11]=k_P1r+k_P1g
	A[11][11]=-1.0*(A[8][11]+A[9][11]+A[10][11]+A[12][11]+A[13][11]+A[15][11])
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

	p=null_space(A)
	total=0.0
	for j in range(16):
		total+=p[j][0]
	for j in range(16):
		p[j][0]=p[j][0]*1.0/total
	p_P1=p[1][0]+p[5][0]+p[6][0]+p[7][0]+p[12][0]+p[13][0]+p[14][0]+p[15][0]
	p_P2=p[2][0]+p[5][0]+p[8][0]+p[9][0]+p[11][0]+p[13][0]+p[14][0]+p[15][0]
	p_N1=p[3][0]+p[6][0]+p[8][0]+p[10][0]+p[11][0]+p[12][0]+p[14][0]+p[15][0]
	p_N2=p[4][0]+p[7][0]+p[9][0]+p[10][0]+p[11][0]+p[12][0]+p[13][0]+p[15][0]

	J1=k_N2l*(1-p_N2)
	J2=k_lN2*p_N2
	J3=k_P1r*(1-p_P1)
	J4=k_rP1*p_P1
	J5=k_P2r*(1-p_P2)
	J6=k_rP2*p_P2
	J7=k_P1g*(1-p_P1)
	J8=k_gP1*p_P1
	J9=k_P2g*(1-p_P2)
	J10=k_gP2*p_P2
	J11=k_N1g*(1-p_N1)
	J12=k_gN1*p_N1

	Jg=J8-J7+J12-J11+J10-J9

	Vout-=1.0*Jg*tint/Cg
	dissipation=1.0*(J1-J2)*(mu_l-mu_g)+1.0*(J3-J4+J5-J6)*(mu_r-mu_g)

	return Vout,dissipation

V_D=5.0
N_NAND=4
V_DATA=0.0
V_WE=0.01
Vin_A_NAND=np.zeros(N_NAND)
Vin_B_NAND=np.zeros(N_NAND)
Vout_NAND=np.zeros(N_NAND)
dissipation_NAND=np.zeros(N_NAND)
Vin_NOT=0.0
Vout_NOT=0.0
dissipation_NOT=0.0
dissipation_tot=0.0

tint=10
T=10000000
Ntot=int(T/tint)

Vin_A_NAND[0]=V_DATA
Vin_B_NAND[0]=V_WE
Vin_A_NAND[1]=V_WE
Vin_NOT = V_DATA

Vout_NAND[2]=0.0
Vout_NAND[3]=0.0

for i in range(Ntot):
    if(i%100==0):
        file=open('result_Dflipflop_D=0_WE=0.01_out=0_0','a')
        file.write("%d	%.10lf   %.10lf	%.10lf	%.10lf	%.10lf	%.10lf	%.10lf	%.10lf	%.10lf	%.10lf\n" %(i*tint,Vout_NOT,Vout_NAND[0],Vout_NAND[1],Vout_NAND[2],Vout_NAND[3],dissipation_NOT,dissipation_NAND[0],dissipation_NAND[1],dissipation_NAND[2],dissipation_NAND[3]))
        file.close()

    Vin_B_NAND[1] = Vout_NOT
    Vin_A_NAND[2] = Vout_NAND[0]
    Vin_B_NAND[2] = Vout_NAND[3]
    Vin_A_NAND[3] = Vout_NAND[2]
    Vin_B_NAND[3] = Vout_NAND[1]

    dissipation_tot = 0.0
    Vout_NOT, dissipation_NOT = NOT_propagation(Vin_NOT,Vout_NOT,V_D)
    dissipation_tot+=dissipation_NOT
    for j in range(N_NAND):
        Vout_NAND[j], dissipation_NAND[j] = NAND_propagation(Vin_A_NAND[j],Vin_B_NAND[j],Vout_NAND[j],V_D)
        dissipation_tot+=dissipation_NAND[j]
