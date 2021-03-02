#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <algorithm>	//std::copy
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

using namespace std;

double fRand(double fMin, double fMax){
	double f = (double)rand()/RAND_MAX;
	return fMin + f*(fMax-fMin);
}

double Fermi(double x){
	return 1.0/(exp(x)+1);
}

double Bose(double x){
	if(x>1e-4){
		return 1.0/(exp(x)-1);
	}
	else{
		//printf("divided by %g in Bose distribution!\n",x);
		return 1e5;
	}
}

double NOT_propagation(int& N_N, int& N_P, double Vin, double& Vout, double& dissipation, double Vd, double tint, gsl_rng* rng)
{
	double Gamma_l=0.2, Gamma_r=0.2, Gamma=0.2, Gamma_g=0.2;
	double E_N, E_P;
	double mu_l=0.0, mu_r=0.0-Vd, mu_g=0.0;
	double kBT=1.0;
	double Cg=200.0;
	double k_Nl,k_lN,k_rP,k_Pr,k_Ng,k_gN,k_Pg,k_gP,k_NP,k_PN;
	double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w,p2;
	double t=0.0,dt;
	int J1,J2,J3,J4;
	int i,j,k;
	double time;

	E_N=1.5*Vd-Vin;
  	E_P=Vin;

	mu_g=0.0-Vout;

	k_Nl=Gamma_l*Fermi((E_N-mu_l)/kBT);
        k_lN=Gamma_l*(1.0-Fermi((E_N-mu_l)/kBT));
        k_rP=Gamma_r*(1.0-Fermi((E_P-mu_r)/kBT));
        k_Pr=Gamma_r*Fermi((E_P-mu_r)/kBT);
        k_Ng=Gamma_g*Fermi((E_N-mu_g)/kBT);
        k_gN=Gamma_g*(1.0-Fermi((E_N-mu_g)/kBT));
        k_Pg=Gamma_g*Fermi((E_P-mu_g)/kBT);
        k_gP=Gamma_g*(1.0-Fermi((E_P-mu_g)/kBT));
        if(E_N>E_P){
		k_NP=Gamma*Bose((E_N-E_P)/kBT);
       		k_PN=Gamma*(1+Bose((E_N-E_P)/kBT));
	}
	else{
                k_PN=Gamma*Bose((E_P-E_N)/kBT);
                k_NP=Gamma*(1+Bose((E_P-E_N)/kBT));
        }

	t=0.0;
	J1=0;
	J2=0;
	J3=0;
	J4=0;

	while(t<tint){
		w1 = k_Nl*(1-N_N);
		w2 = k_lN*N_N;
		w3 = k_rP*N_P;
		w4 = k_Pr*(1-N_P);
                w5 = k_Ng*(1-N_N);
                w6 = k_gN*N_N;
                w7 = k_Pg*(1-N_P);
                w8 = k_gP*N_P;
                w9 = k_PN*N_N*(1-N_P);
                w10 = k_NP*N_P*(1-N_N);
                w = w1+w2+w3+w4+w5+w6+w7+w8+w9+w10;

		dt = -1.0*log(gsl_rng_uniform(rng))/w;
		t = t+dt;

		if(t<tint){
			p2=gsl_rng_uniform(rng);
			if(p2<w1/w){
				N_N=1;
				J3+=1;
			}
			else if(p2<(w1+w2)/w){
				N_N=0;
				J3-=1;
			}
			else if(p2<(w1+w2+w3)/w){
				N_P=0;
				J4+=1;
			}
			else if(p2<(w1+w2+w3+w4)/w){
				N_P=1;
				J4-=1;
			}
			else if(p2<(w1+w2+w3+w4+w5)/w){
				N_N=1;
				J1-=1;
			}
			else if(p2<(w1+w2+w3+w4+w5+w6)/w){
				N_N=0;
				J1+=1;
			}
			else if(p2<(w1+w2+w3+w4+w5+w6+w7)/w){
				N_P=1;
				J2-=1;
			}
			else if(p2<(w1+w2+w3+w4+w5+w6+w7+w8)/w){
				N_P=0;
				J2+=1;
			}
			else if(p2<(w1+w2+w3+w4+w5+w6+w7+w8+w9)/w){
				N_N=0;
				N_P=1;
			}
			else{
				N_N=1;
				N_P=0;
			}
		}
	}

	Vout-=1.0*(J1+J2)/Cg;
	dissipation=1.0*J3*(mu_l-mu_g)-1.0*J4*(mu_r-mu_g);

	return Vout;
	}

double NAND_propagation(int& N_N1, int& N_N2, int& N_P1, int& N_P2, double VinA, double VinB, double& Vout, double& dissipation, double Vd, double tint, gsl_rng* rng)
{
	double Gamma_l=0.2, Gamma_r=0.2, Gamma=0.2, Gamma_g=0.2;
	double mu_l=0.0, mu_r=0.0-Vd, mu_g=0.0;
	double kBT=1.0;
	double Cg=200.0;
	double k_N2l,k_lN2,k_rP1,k_P1r,k_P2r,k_rP2,k_N1g,k_gN1,k_P1g,k_gP1,k_P2g,k_gP2,k_N1N2,k_N2N1,k_N1P1,k_P1N1,k_N1P2,k_P2N1,k_P1P2,k_P2P1;
	int Nreaction=21;
	vector<double> w(Nreaction),wsum(Nreaction);
	double wtot,p2;
	double t=0.0,dt;
	int J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12;
	double E_P1,E_P2,E_N1,E_N2;
	int i,j,k;

	E_P1=VinA;
	E_P2=VinB;
	E_N1=1.5*Vd-VinA;
	E_N2=1.5*Vd-VinB;

	mu_g=0.0-Vout;

	k_N2l=Gamma_l*Fermi((E_N2-mu_l)/kBT);
        k_lN2=Gamma_l*(1.0-Fermi((E_N2-mu_l)/kBT));
        k_rP1=Gamma_r*(1.0-Fermi((E_P1-mu_r)/kBT));
        k_P1r=Gamma_r*Fermi((E_P1-mu_r)/kBT);
	k_P2r=Gamma_r*Fermi((E_P2-mu_r)/kBT);
        k_rP2=Gamma_r*(1.0-Fermi((E_P2-mu_r)/kBT));
        k_N1g=Gamma_g*Fermi((E_N1-mu_g)/kBT);
        k_gN1=Gamma_g*(1.0-Fermi((E_N1-mu_g)/kBT));
        k_P1g=Gamma_g*Fermi((E_P1-mu_g)/kBT);
        k_gP1=Gamma_g*(1.0-Fermi((E_P1-mu_g)/kBT));
	k_P2g=Gamma_g*Fermi((E_P2-mu_g)/kBT);
        k_gP2=Gamma_g*(1.0-Fermi((E_P2-mu_g)/kBT));
	if(E_N1>E_N2){
                k_N1N2=Gamma*Bose((E_N1-E_N2)/kBT);
                k_N2N1=Gamma*(1+Bose((E_N1-E_N2)/kBT));
	}
        else{
                k_N2N1=Gamma*Bose((E_N2-E_N1)/kBT);
                k_N1N2=Gamma*(1+Bose((E_N2-E_N1)/kBT));
	}
	if(E_N1>E_P1){
        	k_N1P1=Gamma*Bose((E_N1-E_P1)/kBT);
        	k_P1N1=Gamma*(1+Bose((E_N1-E_P1)/kBT));
	}
	else{
		k_P1N1=Gamma*Bose((E_P1-E_N1)/kBT);
		k_N1P1=Gamma*(1+Bose((E_P1-E_N1)/kBT));
	}
	if(E_N1>E_P2){
                k_N1P2=Gamma*Bose((E_N1-E_P2)/kBT);
                k_P2N1=Gamma*(1+Bose((E_N1-E_P2)/kBT));
	}
        else{
                k_P2N1=Gamma*Bose((E_P2-E_N1)/kBT);
                k_N1P2=Gamma*(1+Bose((E_P2-E_N1)/kBT));
	}
        if(E_P1>E_P2){
                k_P1P2=Gamma*Bose((E_P1-E_P2)/kBT);
                k_P2P1=Gamma*(1+Bose((E_P1-E_P2)/kBT));
	}
        else{
                k_P2P1=Gamma*Bose((E_P2-E_P1)/kBT);
                k_P1P2=Gamma*(1+Bose((E_P2-E_P1)/kBT));
	}

	t=0.0;
	J1=0;
	J2=0;
	J3=0;
	J4=0;
	J5=0;
	J6=0;
	J7=0;
	J8=0;
	J9=0;
	J10=0;
	J11=0;
	J12=0;

	while(t<tint){
		w[0]=0.0;
		w[1] = k_N2l*(1-N_N2);
		w[2] = k_lN2*N_N2;
		w[3] = k_P1r*(1-N_P1);
		w[4] = k_rP1*N_P1;
		w[5] = k_P2r*(1-N_P2);
		w[6] = k_rP2*N_P2;
		w[7] = k_P1g*(1-N_P1);
		w[8] = k_gP1*N_P1;
		w[9] = k_P2g*(1-N_P2);
		w[10] = k_gP2*N_P2;
    w[11] = k_N1g*(1-N_N1);
    w[12] = k_gN1*N_N1;
    w[13]=k_N1N2*N_N2*(1-N_N1);
    w[14]=k_N2N1*N_N1*(1-N_N2);
    w[15]=k_N1P1*N_P1*(1-N_N1);
    w[16]=k_P1N1*N_N1*(1-N_P1);
    w[17]=k_N1P2*N_P2*(1-N_N1);
    w[18]=k_P2N1*N_N1*(1-N_P2);
    w[19]=k_P1P2*N_P2*(1-N_P1);
    w[20]=k_P2P1*N_P1*(1-N_P2);

		wsum[0]=w[0];
		for(j=1;j<Nreaction;j++){
				wsum[j]=wsum[j-1]+w[j];
		}
    wtot=wsum[Nreaction-1];

		dt = -1.0*log(gsl_rng_uniform(rng))/wtot;
		t = t+dt;

		if(t<tint){
			p2=gsl_rng_uniform(rng);
			if(p2<wsum[1]/wtot){
				N_N2=1;
				J1+=1;
			}
			else if(p2<wsum[2]/wtot){
				N_N2=0;
				J2+=1;
			}
			else if(p2<wsum[3]/wtot){
				N_P1=1;
				J3+=1;
			}
			else if(p2<wsum[4]/wtot){
				N_P1=0;
				J4+=1;
			}
			else if(p2<wsum[5]/wtot){
				N_P2=1;
				J5+=1;
			}
			else if(p2<wsum[6]/wtot){
				N_P2=0;
				J6+=1;
			}
			else if(p2<wsum[7]/wtot){
				N_P1=1;
				J7+=1;
			}
			else if(p2<wsum[8]/wtot){
				N_P1=0;
				J8+=1;
			}
			else if(p2<wsum[9]/wtot){
				N_P2=1;
				J9+=1;
			}
			else if(p2<wsum[10]/wtot){
				N_P2=0;
				J10+=1;
			}
			else if(p2<wsum[11]/wtot){
				N_N1=1;
				J11+=1;
			}
			else if(p2<wsum[12]/wtot){
				N_N1=0;
				J12+=1;
			}
			else if(p2<wsum[13]/wtot){
				N_N1=1;
				N_N2=0;
			}
			else if(p2<wsum[14]/wtot){
				N_N1=0;
				N_N2=1;
			}
			else if(p2<wsum[15]/wtot){
				N_N1=1;
				N_P1=0;
			}
			else if(p2<wsum[16]/wtot){
				N_N1=0;
				N_P1=1;
			}
			else if(p2<wsum[17]/wtot){
				N_N1=1;
				N_P2=0;
			}
			else if(p2<wsum[18]/wtot){
				N_N1=0;
				N_P2=1;
			}
			else if(p2<wsum[19]/wtot){
				N_P1=1;
				N_P2=0;
			}
			else{
				N_P1=0;
				N_P2=1;
			}
		}
	}
	Vout-=1.0*(J8-J7+J12-J11+J10-J9)/Cg;
	dissipation=1.0*(J1-J2)*(mu_l-mu_g)+1.0*(J3-J4+J5-J6)*(mu_r-mu_g);

  return Vout;
	}

int main()
{
  srand(time(0));
	//initializing the MPI environment
	MPI_Init(NULL,NULL);
	int nprocs,rank;
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status status;

	gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL)+rank);

	int i,j,k,s,index;

	int Nw = 25600;
	int n = Nw/nprocs;
	double tint = 10.0;
	double Tcycle = 1000000.0;
	int Ntint = 100000;
	int Ncycle = 200;

	int Ngate = 4;
	double Vd = 8.0;
	vector< vector<int> > N_N1(n, vector<int>(Ngate));
	vector< vector<int> > N_N2(n, vector<int>(Ngate));
	vector< vector<int> > N_P1(n, vector<int>(Ngate));
	vector< vector<int> > N_P2(n, vector<int>(Ngate));
	vector<int> N_N(n);
	vector<int> N_P(n);
	vector< vector<double> > Vin_A(n, vector<double>(Ngate));
	vector< vector<double> > Vin_B(n, vector<double>(Ngate));
	vector< vector<double> > Vout(n, vector<double>(Ngate));
	vector<double> Vin_NOT(n);
	vector<double> Vout_NOT(n);
	vector<double> output_local_NOT(n);
	vector<double> output_global_NOT(Nw);
	vector<double> dissipation_local_NOT(n);
	vector< vector<double> > output_local(Ngate, vector<double>(n));
	vector< vector<double> > output_global(Ngate, vector<double>(Nw));
	vector< vector<double> > dissipation_local(Ngate, vector<double>(n));
	vector<double> dissipation_tot_local(n);
	vector<double> dissipation_tot_global(Nw);
	vector<double> output_avg(Ncycle);
	vector<int> error_count(Ncycle);
	vector<double> dissipation_avg(Ncycle);

	FILE *result;
	char path[100];

	int Nbin_V=1800;
  double Vmin=-0.49;
  double Vmax=8.51;
  vector<double> bin_V(Nbin_V);
  vector<int> hist_V(Nbin_V);

  for(i=0;i<Nbin_V;i++){
      bin_V[i]=Vmin+i*0.005;
  }

	for(j=0;j<n;j++){
		Vout_NOT[j] =  0.0002;
		Vout[j][0] = 0.0002;
		Vout[j][1] = 0.0002;
		Vout[j][2] = 0.0002;
		Vout[j][3] = 0.0002;
		N_N[j] = 0;
		N_P[j] = 0;
		for(k=0;k<Ngate;k++){
			N_N1[j][k] = 0;
			N_N2[j][k] = 0;
			N_P1[j][k] = 0;
			N_P2[j][k] = 0;
		}
	}

	for(s=0;s<Ncycle;s++){
		// first half cycle
		for(j=0;j<n;j++){
			if(s%2==0){
				Vin_A[j][0] = Vd;
				Vin_NOT[j] = Vd;
			}
			else{
				Vin_A[j][0] = 0.0;
				Vin_NOT[j] = 0.0;
			}
			Vin_B[j][0] = Vd+0.01;
			Vin_A[j][1] = Vd+0.01;
			dissipation_tot_local[j] = 0.0;
		}

		for(i=0;i<Ntint/2;i++){
		for(j=0;j<n;j++){
			Vin_B[j][1] = Vout_NOT[j];
			Vin_A[j][2] = Vout[j][0];
			Vin_B[j][2] = Vout[j][3];
			Vin_A[j][3] = Vout[j][2];
			Vin_B[j][3] = Vout[j][1];
			output_local_NOT[j] = NOT_propagation(N_N[j],N_P[j],Vin_NOT[j],Vout_NOT[j],dissipation_local_NOT[j],Vd,tint,rng);
			dissipation_tot_local[j]+=dissipation_local_NOT[j];
			for(k=0;k<Ngate;k++){
			output_local[k][j] = NAND_propagation(N_N1[j][k],N_N2[j][k],N_P1[j][k],N_P2[j][k],Vin_A[j][k],Vin_B[j][k],Vout[j][k],dissipation_local[k][j],Vd,tint,rng);
			dissipation_tot_local[j]+=dissipation_local[k][j];
			}
		}
		}

		// second half Ncycle
		for(j=0;j<n;j++){
			Vin_B[j][0] = 0.01;
			Vin_A[j][1] = 0.01;
		}

		for(i=0;i<Ntint/2;i++){
		for(j=0;j<n;j++){
			Vin_B[j][1] = Vout_NOT[j];
			Vin_A[j][2] = Vout[j][0];
			Vin_B[j][2] = Vout[j][3];
			Vin_A[j][3] = Vout[j][2];
			Vin_B[j][3] = Vout[j][1];
			output_local_NOT[j] = NOT_propagation(N_N[j],N_P[j],Vin_NOT[j],Vout_NOT[j],dissipation_local_NOT[j],Vd,tint,rng);
			dissipation_tot_local[j]+=dissipation_local_NOT[j];
			for(k=0;k<Ngate;k++){
			output_local[k][j] = NAND_propagation(N_N1[j][k],N_N2[j][k],N_P1[j][k],N_P2[j][k],Vin_A[j][k],Vin_B[j][k],Vout[j][k],dissipation_local[k][j],Vd,tint,rng);
			dissipation_tot_local[j]+=dissipation_local[k][j];
			}
		}
		}

		MPI_Gather(&output_local[2][0],n,MPI_DOUBLE,&output_global[2][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&dissipation_tot_local[0],n,MPI_DOUBLE,&dissipation_tot_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
		if(rank==0){
			error_count[s]=0;
			output_avg[s]=0.0;
			dissipation_avg[s]=0.0;
			for(j=0;j<Nw;j++){
				output_avg[s]+=output_global[2][j];
				dissipation_avg[s]+=dissipation_tot_global[j];
				if(s%2==0 && output_global[2][j]<(0.98*Vd-1e-4)){
					error_count[s]+=1;
				}
				if(s%2==1 && output_global[2][j]>(0.02*Vd+1e-4)){
					error_count[s]+=1;
				}
			}

			result=fopen("Vavg","a");
			fprintf(result,"%d	%lf	%lf	%lf\n",s,output_avg[s]/Nw,1.0*error_count[s]/Nw,dissipation_avg[s]/Nw);
			fclose(result);

			if((s+1)%10==0){
			for(j=0;j<Nbin_V-1;j++){
          hist_V[j]=0;
      }
      for(k=0;k<Nw;k++){
          for(j=0;j<Nbin_V-1;j++){
              if(output_global[2][k]>=bin_V[j] && output_global[2][k]<bin_V[j+1]){
                  hist_V[j]+=1;
              }
          }
      }
      sprintf(path,"hist_V_%d",s+1);
      result=fopen(path,"w");
      for(j=0;j<Nbin_V;j++){
          fprintf(result,"%lf     %d\n",bin_V[j],hist_V[j]);
      }
      fclose(result);
			}

			if(s>=50){
			result=fopen("/global/scratch/chloegao/Vtraj_Vd=8_in=10_T=1e6","a");
				for(j=0;j<Nw;j++){
					fprintf(result,"%lf	",output_global[2][j]);
				}
				fprintf(result,"\n");
				fclose(result);
			}
		}
	}

	MPI_Finalize();
	return 0;
}
