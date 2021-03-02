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
		printf("divided by %g in Bose distribution!\n",x);
		return 1e5;
	}
}

double NAND_propagation(int& N_N1, int& N_N2, int& N_P1, int& N_P2, double VinA, double VinB, double& Vout, double Vd, double& dissipation, double tint, gsl_rng* rng)
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

	int Nw = 2400000;
	int n = Nw/nprocs;
	double tint = 10.0, tobs = 6000000.0;
	int Ntint = tobs/tint;
	int Ngate = 2;
	double Vd = 5.0;
	vector< vector<int> > N_N1(n, vector<int>(Ngate));
	vector< vector<int> > N_N2(n, vector<int>(Ngate));
	vector< vector<int> > N_P1(n, vector<int>(Ngate));
	vector< vector<int> > N_P2(n, vector<int>(Ngate));
	vector< vector<double> > Vin_A(n, vector<double>(Ngate));
	vector< vector<double> > Vin_B(n, vector<double>(Ngate));
	vector< vector<double> > Vout(n, vector<double>(Ngate));
	vector< vector<double> > output_local(Ngate, vector<double>(n));
	vector< vector<double> > output_global(Ngate, vector<double>(Nw));
	vector< vector<double> > dissipation_local(Ngate, vector<double>(n));
	vector< vector<double> > dissipation_global(Ngate, vector<double>(Nw));
	vector<double> output_avg(Ngate), dissipation_avg(Ngate);
	double t_avg;

	FILE *result;
	char path[100];

	int Nbin_V=1200;
	double V1min=-0.49;
	double V1max=5.51;
	vector<double> bin_V1(Nbin_V);
	vector<int> hist_V1(Nbin_V);

	double V2min=-0.49;
	double V2max=5.51;
	vector<double> bin_V2(Nbin_V);
	vector<int> hist_V2(Nbin_V);


	for(j=0;j<Nbin_V;j++){
		bin_V1[j]=V1min+j*0.005;
		bin_V2[j]=V2min+j*0.005;
	}

	for(j=0;j<n;j++){
		Vin_A[j][0] = Vd;
		Vin_B[j][1] = Vd;
		Vout[j][0] = 5.0002;
		Vout[j][1] = 0.0002;
		for(k=0;k<Ngate;k++){
			N_N1[j][k] = 0;
			N_N2[j][k] = 0;
			N_P1[j][k] = 0;
			N_P2[j][k] = 0;
		}
	}

	for(i=0;i<Ntint;i++){
		for(j=0;j<n;j++){
			Vin_A[j][1] = Vout[j][0];
			Vin_B[j][0] = Vout[j][1];
			for(k=0;k<Ngate;k++){
				output_local[k][j] = NAND_propagation(N_N1[j][k],N_N2[j][k],N_P1[j][k],N_P2[j][k],Vin_A[j][k],Vin_B[j][k],Vout[j][k],Vd,dissipation_local[k][j],tint,rng);
			}
		}
		if((i+1)%10000==0){
		  for(k=0;k<Ngate;k++){
	    	MPI_Gather(&output_local[k][0],n,MPI_DOUBLE,&output_global[k][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
      	MPI_Gather(&dissipation_local[k][0],n,MPI_DOUBLE,&dissipation_global[k][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    	}
			if(rank==0){
				for(k=0;k<Ngate;k++){
					output_avg[k] = 0.0;
					dissipation_avg[k] = 0.0;
				}
				for(j=0;j<Nw;j++){
					for(k=0;k<Ngate;k++){
						output_avg[k] += output_global[k][j]/Nw;
						dissipation_avg[k] += dissipation_global[k][j]/Nw;
					}
				}
				result=fopen("time_result","a");
				fprintf(result,"%lf	%lf	%lf	%lf	%lf\n",(i+1)*tint,output_avg[0],output_avg[1],dissipation_avg[0],dissipation_avg[1]);
				fclose(result);

				for(j=0;j<Nbin_V-1;j++){
                                        hist_V1[j]=0;
                                        hist_V2[j]=0;
                                }

				for(k=0;k<Nw;k++){
					for(j=0;j<Nbin_V-1;j++){
						if(output_global[0][k]>=bin_V1[j] && output_global[0][k]<bin_V1[j+1]){
							hist_V1[j]+=1;
						}
						if(output_global[1][k]>=bin_V2[j] && output_global[1][k]<bin_V2[j+1]){
							hist_V2[j]+=1;
						}
					}
				}

				sprintf(path,"hist_V1_%d",i+1);
				result=fopen(path,"w");
		    for(j=0;j<Nbin_V;j++){
		      fprintf(result,"%lf     %d\n",bin_V1[j],hist_V1[j]);
		    }
		    fclose(result);

				sprintf(path,"hist_V2_%d",i+1);
				result=fopen(path,"w");
		    for(j=0;j<Nbin_V;j++){
		      fprintf(result,"%lf     %d\n",bin_V2[j],hist_V2[j]);
		    }
		    fclose(result);
			}
		}
	}


	if(rank==0){
		printf("average Vg is %lf,%lf, average delay is %lf, average dissipation rate is %lf,%lf\n",output_avg[0],output_avg[1],t_avg,dissipation_avg[0],dissipation_avg[1]);
	}

	MPI_Finalize();
	return 0;
}
