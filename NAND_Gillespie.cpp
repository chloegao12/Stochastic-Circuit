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

struct Walker_info{
	vector<double> x;
	vector<double> v;
	int parent;
	double Qavg;
};

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
	else	return 1e5;
}

// this function takes the initial value x&v and modify it
// and returns the cumulative current Q within tint
double propagation(double& delay, double& dissipation_steady, double& dissipation_nonsteady, double tint, int Nint, gsl_rng* rng)
{
	double Gamma_l=0.2, Gamma_r=0.2, Gamma=0.2, Gamma_g=0.2;
	double E_N,E_P;
	double Vd=5.0;
	double VinA=5.0, VinB=0.0;
	double mu_l=0.0, mu_r=0.0-Vd, mu_g=0.0;
	double Vout=0.0;
	double kBT=1.0;
	double Cg=200.0;
	double k_N2l,k_lN2,k_rP1,k_P1r,k_P2r,k_rP2,k_N1g,k_gN1,k_P1g,k_gP1,k_P2g,k_gP2,k_N1N2,k_N2N1,k_N1P1,k_P1N1,k_N1P2,k_P2N1,k_P1P2,k_P2P1;
	//double w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,
	int Nreaction=21;
	vector<double> w(Nreaction),wsum(Nreaction);
	double wtot,p2;
	double t=0.0,dt;
	int J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12;
	int N_N1=0,N_N2=0,N_P1=0,N_P2=0;
	double E_P1,E_P2,E_N1,E_N2;
	int i,j,k;
	vector<double> dissipation(Nint),output(Nint);
	double time;

	E_P1=VinA;
	E_P2=VinB;
	E_N1=1.5*Vd-VinA;
	E_N2=1.5*Vd-VinB;

	for(i=0;i<Nint;i++){
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
	output[i]=Vout;
	dissipation[i]=1.0*(J1-J2)*(mu_l-mu_g)+1.0*(J3-J4+J5-J6)*(mu_r-mu_g);
	}

	for(i=0;i<Nint;i++){
		if(output[i]>(0.98*Vd-1e-4)){
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
