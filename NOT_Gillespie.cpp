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
	return 1.0/(exp(x)-1);
}

// this function takes the initial value x&v and modify it
// and returns the cumulative current Q within tint
double propagation(double& delay, double& dissipation_steady, double& dissipation_nonsteady, double tint, int Nint, gsl_rng* rng)
{
	double Gamma_l=0.2, Gamma_r=0.2, Gamma=0.2, Gamma_g=0.2;
	double E_N,E_P;
	double Vd=5.0;
	double Vin;
	double mu_l=0.0, mu_r=0.0-Vd, mu_g=0.0;
	double Vout=0.0;
	double kBT=1.0;
	double Cg=200.0;
	double k_Nl,k_lN,k_rP,k_Pr,k_Ng,k_gN,k_Pg,k_gP,k_NP,k_PN;
	double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w,p2;
	double t=0.0,dt;
	int J1,J2,J3,J4;
	int N_N=0,N_P=0;
	int i,j,k;
	vector<double> dissipation(Nint),output(Nint);
	double time;

	E_N=1.5*Vd-Vin;
	E_P=Vin;

	for(i=0;i<Nint;i++){
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
	output[i]=Vout;
	dissipation[i]=1.0*J3*(mu_l-mu_g)-1.0*J4*(mu_r-mu_g);
	}

	for(i=0;i<Nint;i++){
		if(output[i]>=0.98*Vd){
			time=i;
			break;
		}
	}
	delay=time*tint;
	dissipation_nonsteady=0.0;
	for(i=0;i<(time+1);i++){
		dissipation_nonsteady+=dissipation[i];
	}
	dissipation_steady=0.0;
	for(i=500000;i<Nint;i++){
		dissipation_steady+=dissipation[i];
	}
	dissipation_steady=dissipation_steady/100000/tint;

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

	//parameters for CANSS
	int Nw = 24000000;
	int n = Nw/nprocs;
	double tint = 10.0, tobs = 6000000.0;
	int Ntint = tobs/tint;
	vector<double> output_local(n),output_global(Nw);
	vector<double> t_local(n),t_global(Nw);
	vector<double> dissipation_nonsteady_local(n),dissipation_nonsteady_global(Nw);
	vector<double> dissipation_steady_local(n),dissipation_steady_global(Nw);
	double output_avg,t_avg,dissipation_nonsteady_avg,dissipation_steady_avg;

	FILE *result;

	for(j=0;j<n;j++){
		output_local[j]=propagation(t_local[j],dissipation_steady_local[j],dissipation_nonsteady_local[j],tint,Ntint,rng);
	}

	MPI_Gather(&t_local[0],n,MPI_DOUBLE,&t_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&output_local[0],n,MPI_DOUBLE,&output_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&dissipation_steady_local[0],n,MPI_DOUBLE,&dissipation_steady_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&dissipation_nonsteady_local[0],n,MPI_DOUBLE,&dissipation_nonsteady_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(rank==0){
		output_avg=0.0;
		t_avg=0.0;
		dissipation_nonsteady_avg=0.0;
		dissipation_steady_avg=0.0;

		for(j=0;j<Nw;j++){
			output_avg+=output_global[j];
			t_avg+=t_global[j];
			dissipation_nonsteady_avg+=dissipation_nonsteady_global[j];
			dissipation_steady_avg+=dissipation_steady_global[j];
		}
		output_avg=output_avg/Nw;
		t_avg=t_avg/Nw;
		dissipation_nonsteady_avg=dissipation_nonsteady_avg/Nw;
		dissipation_steady_avg=dissipation_steady_avg/Nw;

		printf("average Vg is %lf, average delay is %lf, average transient dissipation is %lf, average dissipation rate in steady state is %lf\n",output_avg,t_avg,dissipation_nonsteady_avg,dissipation_steady_avg);


		int Nbin_t=120000;
		double tmin=2.0;
		double tmax=4800002.0;
		vector<double> bin_t(Nbin_t);
		vector<int> hist_t(Nbin_t);

		int Nbin_V=300;
		double Vmin=4.01;
		double Vmax=5.51;
		vector<double> bin_V(Nbin_V);
		vector<int> hist_V(Nbin_V);

		int Nbin_dis_nonsteady=600;
        	double dmin_nonsteady=2001;
        	double dmax_nonsteady=5001;
        	vector<double> bin_dis_nonsteady(Nbin_dis_nonsteady);
        	vector<int> hist_dis_nonsteady(Nbin_dis_nonsteady);

		int Nbin_dis_steady=60;
        	double dmin_steady=0.0001;
        	double dmax_steady=0.0007;
        	vector<double> bin_dis_steady(Nbin_dis_steady);
        	vector<int> hist_dis_steady(Nbin_dis_steady);

		for(i=0;i<Nbin_t;i++){
			bin_t[i]=tmin+i*40;
		}

		for(i=0;i<Nbin_V;i++){
                	bin_V[i]=Vmin+i*0.005;
        	}

		for(i=0;i<Nbin_dis_nonsteady;i++){
                	bin_dis_nonsteady[i]=dmin_nonsteady+i*5;
        	}

		for(i=0;i<Nbin_dis_steady;i++){
                	bin_dis_steady[i]=dmin_steady+i*0.00001;
        	}

		for(k=0;k<Nw;k++){
			for(j=0;j<Nbin_t-1;j++){
				if(t_global[k]>=bin_t[j] && t_global[k]<bin_t[j+1]){
					hist_t[j]+=1;
				}
			}

			for(j=0;j<Nbin_V-1;j++){
				if(output_global[k]>=bin_V[j] && output_global[k]<bin_V[j+1]){
					hist_V[j]+=1;
				}
			}

			for(j=0;j<Nbin_dis_nonsteady-1;j++){
          			if(dissipation_nonsteady_global[k]>=bin_dis_nonsteady[j] && dissipation_nonsteady_global[k]<bin_dis_nonsteady[j+1]){
                        		hist_dis_nonsteady[j]+=1;
                		}
        		}

			for(j=0;j<Nbin_dis_steady-1;j++){
         	        	if(dissipation_steady_global[k]>=bin_dis_steady[j] && dissipation_steady_global[k]<bin_dis_steady[j+1]){
                        		hist_dis_steady[j]+=1;
                		}
        		}
		}

		result=fopen("hist_t_tint=10_C=200","w");
		for(i=0;i<Nbin_t;i++){
		fprintf(result,"%lf	%d\n",bin_t[i],hist_t[i]);
		}
		fclose(result);

		result=fopen("hist_V_tint=10_C=200","w");
        	for(i=0;i<Nbin_V;i++){
                fprintf(result,"%lf     %d\n",bin_V[i],hist_V[i]);
        	}
        	fclose(result);

		result=fopen("hist_dis_nonsteady_tint=10_C=200","w");
        	for(i=0;i<Nbin_dis_nonsteady;i++){
                fprintf(result,"%lf     %d\n",bin_dis_nonsteady[i],hist_dis_nonsteady[i]);
        	}
        	fclose(result);

		result=fopen("hist_dis_steady_tint=10_C=200","w");
        	for(i=0;i<Nbin_dis_steady;i++){
                fprintf(result,"%lf     %d\n",bin_dis_steady[i],hist_dis_steady[i]);
        	}
        	fclose(result);

	}

	MPI_Finalize();
	return 0;
}
