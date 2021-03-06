from pylab import *
import numpy as np
import scipy
from scipy.linalg import null_space

def Fermi(x):
        return 1.0/(exp(x)+1)

def Bose(x):
        return 1.0/(exp(x)-1)

A=zeros((4,4))
c=zeros(4)

Gamma_l=0.2
Gamma_r=0.2
Gamma=0.2
Gamma_g=0.2

V_D=5.0
mu_l=0.0
mu_r=mu_l-V_D
kBT=1.0
Cg=200.0

Ngate=40
Vin=np.zeros(Ngate)
Vout=np.zeros(Ngate)
current=np.zeros(Ngate)
dissipation=np.zeros(Ngate)

Vin[0]=0.0

tint=10
T=10000000
Ntot=int(T/tint)

for i in range(Ntot):
	for k in range(Ngate):
		if(k>0):Vin[k]=Vout[k-1]
	for k in range(Ngate):
		E_N=1.5*V_D-Vin[k]
		E_P=Vin[k]

		mu_g=0.0-Vout[k]

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

		Vout[k]-=1.0*(J1+J2)*tint/Cg
		current[k]=(J1+J2)*tint
		dissipation[k]=1.0*J3*tint*(mu_l-mu_g)-1.0*J4*tint*(mu_r-mu_g)
	
	if(i%10==0):
		file=open('result_Vd=5_N=40','a')	
		file.write("%d	%g	" %(i*tint,current[Ngate-1]))
		for k in range(Ngate):
			file.write("%.7lf	" %Vout[k])
		file.write("\n")
		file.close()

#file=open('result_EN=0.692911_EP=4.306989_Vd=5','w')
#for i in range(Ntot):
#        file.write("%lf %lf     %lf      %g\n" %(i*tint,output[i],current[i],dissipation[i]))
#file.close()
