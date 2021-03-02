from pylab import *
import numpy as np
import scipy
from scipy import integrate

def Gauss_in0(x):
	sigma2 = 0.005 
	mu = 5.0159655
	return 1/sqrt(2*np.pi*sigma2)*exp(-1.0/2/sigma2*(x-mu)*(x-mu))

def Gauss_in1(x):
	sigma2 = 0.005
	mu = 0.090022
	return 1/sqrt(2*np.pi*sigma2)*exp(-1.0/2/sigma2*(x-mu)*(x-mu))

Vd=5.0
p_cond=np.zeros(6)
p_joint=np.zeros(6)
p_X=np.zeros(2)
p_Y=np.zeros(3)

p_cond[0]=integrate.quad(Gauss_in0,-1*np.inf,0.02*Vd,epsabs=1e-28)[0]
p_cond[1]=integrate.quad(Gauss_in0,0.98*Vd,np.inf,epsabs=1e-28)[0]
p_cond[2]=1.0-p_cond[0]-p_cond[1]
p_cond[3]=integrate.quad(Gauss_in1,-1*np.inf,0.02*Vd,epsabs=1e-28)[0]
p_cond[4]=integrate.quad(Gauss_in1,0.98*Vd,np.inf,epsabs=1e-28)[0]
p_cond[5]=1.0-p_cond[3]-p_cond[4]
print(p_cond)

pis=np.arange(0.0001,1.0,0.0001)
Cs=[]

for pi in pis:
                p_X[0]=pi
                p_X[1]=1.0-pi
                for i in range(3):
                        p_joint[i]=p_cond[i]*p_X[0]
                for i in range(3,6):
                        p_joint[i]=p_cond[i]*p_X[1]
                p_Y[0]=p_joint[0]+p_joint[3]
                p_Y[1]=p_joint[1]+p_joint[4]
                p_Y[2]=p_joint[2]+p_joint[5]
                I_mutual=0.0
                if(p_joint[0]>0):
                        I_mutual+=p_joint[0]*np.log2(p_cond[0]/p_Y[0])
                I_mutual+=p_joint[1]*np.log2(p_cond[1]/p_Y[1])
                I_mutual+=p_joint[2]*np.log2(p_cond[2]/p_Y[2])
                I_mutual+=p_joint[3]*np.log2(p_cond[3]/p_Y[0])
                if(p_joint[4]>0):
                        I_mutual+=p_joint[4]*np.log2(p_cond[4]/p_Y[1])
                if(p_joint[5]>0):
                        I_mutual+=p_joint[5]*np.log2(p_cond[5]/p_Y[2])
                Cs=np.append(Cs,I_mutual)
                #print(pi,I_mutual)

print(np.amax(Cs))
