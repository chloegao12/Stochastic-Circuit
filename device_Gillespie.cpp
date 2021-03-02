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

double Fermi(double x){
	return 1.0/(exp(x)+1);
}

double Bose(double x){
	if(x>1e-3)	return 1.0/(exp(x)-1);
	else	return 1e4;
}

void NOT_propagation(int& N_N, int& N_P, double Vin, double& Vout, double& dissipation, double Vd, double tint, gsl_rng* rng)
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
	return;
	}

void NAND_propagation(int& N_N1, int& N_N2, int& N_P1, int& N_P2, double VinA, double VinB, double& Vout, double& dissipation, double Vd, double tint, gsl_rng* rng)
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
	return;
}

double XOR_propagation(vector<int>& N_N1, vector<int>& N_N2, vector<int>& N_P1, vector<int>& N_P2, vector<double>& Vout, vector<double>& NAND_dissipation, double& dissipation_tot, double VA, double VB, double Vd, double tint, gsl_rng* rng)
{
	vector<double> VinA(4),VinB(4);
	VinA[0]=VA;
	VinB[0]=VB;
	VinA[1]=VA;
	VinB[2]=VB;

	VinB[1]=Vout[0];
	VinA[2]=Vout[0];
	VinA[3]=Vout[1];
	VinB[3]=Vout[2];
	for(int j=0;j<4;j++){
		NAND_propagation(N_N1[j],N_N2[j],N_P1[j],N_P2[j],VinA[j],VinB[j],Vout[j],NAND_dissipation[j],Vd,tint,rng);
		dissipation_tot+=NAND_dissipation[j];
	}

	return Vout[3];
}

double Dflip_propagation(vector<int>& N_N1, vector<int>& N_N2, vector<int>& N_P1, vector<int>& N_P2, vector<double>& Vout_NAND, vector<double>& dissipation_NAND, int& N_N, int& N_P, double& Vout_NOT, double& dissipation_NOT, double& dissipation_tot, double VD, double WE, double Vd, double tint, gsl_rng* rng)
{
	int j;
	vector<double> VinA(4),VinB(4);
	double Vin;
	VinA[0]=VD;
	VinB[0]=WE;
	VinA[1]=WE;
	Vin=VD;

	VinB[1]=Vout_NOT;
	VinA[2]=Vout_NAND[0];
	VinB[2]=Vout_NAND[3];
	VinA[3]=Vout_NAND[2];
	VinB[3]=Vout_NAND[1];
	NOT_propagation(N_N,N_P,Vin,Vout_NOT,dissipation_NOT,Vd,tint,rng);
	dissipation_tot+=dissipation_NOT;
	for(j=0;j<4;j++){
		NAND_propagation(N_N1[j],N_N2[j],N_P1[j],N_P2[j],VinA[j],VinB[j],Vout_NAND[j],dissipation_NAND[j],Vd,tint,rng);
		dissipation_tot+=dissipation_NAND[j];
	}

	return Vout_NAND[2];
}

int XOR_gate(int x1, int x2)
{
	if((x1+x2)%2==0)	return 0;
	else	return 1;
}

double undigitalize(int x, double Vd)
{
	if(x==0)	return 0.0;
	else	return Vd;
}

int digitalize(double x, double Vd)
{
	if(x>(0.98*Vd-1e-4))	return 1;
	else if(x<(0.02*Vd+1e-4))	return 0;
	else return -1;
}

int main()
{
  srand(time(0));
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

	int Nw = 12800;
	int n = Nw/nprocs;
	double tint = 10.0;
	double Tcycle = 10000000.0;
	int Ntint = Tcycle/tint;

	double Vd = 8.0;
	int NXOR = 2;
	int NDflip = 4;

	// data for 2 XOR gates
	vector< vector <vector<int> > > XOR_N_N1,XOR_N_N2,XOR_N_P1,XOR_N_P2;
	vector< vector <vector<double> > > XOR_Vout,XOR_NAND_dissipation;
	XOR_N_N1.resize(n);
	XOR_N_N2.resize(n);
	XOR_N_P1.resize(n);
	XOR_N_P2.resize(n);
	XOR_Vout.resize(n);
	XOR_NAND_dissipation.resize(n);
	for(i=0;i<n;i++){
		XOR_N_N1.at(i).resize(NXOR);
		XOR_N_N2.at(i).resize(NXOR);
		XOR_N_P1.at(i).resize(NXOR);
		XOR_N_P2.at(i).resize(NXOR);
		XOR_Vout.at(i).resize(NXOR);
		XOR_NAND_dissipation.at(i).resize(NXOR);
	}
	for(i=0;i<n;i++){
		for(j=0;j<NXOR;j++){
			XOR_N_N1.at(i).at(j).resize(4);
			XOR_N_N2.at(i).at(j).resize(4);
			XOR_N_P1.at(i).at(j).resize(4);
			XOR_N_P2.at(i).at(j).resize(4);
			XOR_Vout.at(i).at(j).resize(4);
			XOR_NAND_dissipation.at(i).at(j).resize(4);
		}
	}
	vector< vector<double> > XOR_VA(n,vector<double>(NXOR));
	vector< vector<double> > XOR_VB(n,vector<double>(NXOR));
	vector< vector<double> > XOR_dissipation(n,vector<double>(NXOR));
	vector< vector<double> > output_XOR_local(n, vector<double>(NXOR));

	// data for 4 Dflipflops
	vector< vector <vector<int> > > Dflip_N_N1,Dflip_N_N2,Dflip_N_P1,Dflip_N_P2;
	vector< vector <vector<double> > > Dflip_NAND_Vout,Dflip_NAND_dissipation;
	vector< vector<int> > Dflip_N_N,Dflip_N_P;
	vector< vector<double> > Dflip_NOT_Vout,Dflip_NOT_dissipation;
	Dflip_N_N1.resize(n);
	Dflip_N_N2.resize(n);
	Dflip_N_P1.resize(n);
	Dflip_N_P2.resize(n);
	Dflip_NAND_Vout.resize(n);
	Dflip_NAND_dissipation.resize(n);
	Dflip_N_N.resize(n);
	Dflip_N_P.resize(n);
	Dflip_NOT_Vout.resize(n);
	Dflip_NOT_dissipation.resize(n);
	for(i=0;i<n;i++){
		Dflip_N_N1.at(i).resize(NDflip);
		Dflip_N_N2.at(i).resize(NDflip);
		Dflip_N_P1.at(i).resize(NDflip);
		Dflip_N_P2.at(i).resize(NDflip);
		Dflip_NAND_Vout.at(i).resize(NDflip);
		Dflip_NAND_dissipation.at(i).resize(NDflip);
		Dflip_N_N.at(i).resize(NDflip);
		Dflip_N_P.at(i).resize(NDflip);
		Dflip_NOT_Vout.at(i).resize(NDflip);
		Dflip_NOT_dissipation.at(i).resize(NDflip);
	}
	for(i=0;i<n;i++){
		for(j=0;j<NDflip;j++){
			Dflip_N_N1.at(i).at(j).resize(4);
			Dflip_N_N2.at(i).at(j).resize(4);
			Dflip_N_P1.at(i).at(j).resize(4);
			Dflip_N_P2.at(i).at(j).resize(4);
			Dflip_NAND_Vout.at(i).at(j).resize(4);
			Dflip_NAND_dissipation.at(i).at(j).resize(4);
		}
	}
	vector< vector<double> > Dflip_VD(n,vector<double>(NDflip));
	vector< vector<double> > Dflip_WE(n,vector<double>(NDflip));
	vector< vector<double> > Dflip_dissipation(n,vector<double>(NDflip));
	vector< vector<double> > output_Dflip_local(n, vector<double>(NDflip));

	FILE *result;

	// initialization from zero charges
	for(i=0;i<n;i++){
		for(j=0;j<NXOR;j++){
				XOR_dissipation[i][j]=0.0;
				for(k=0;k<4;k++){
					XOR_N_N1[i][j][k]=0;
					XOR_N_N2[i][j][k]=0;
					XOR_N_P1[i][j][k]=0;
					XOR_N_P2[i][j][k]=0;
					XOR_Vout[i][j][k]=0.0;
					XOR_NAND_dissipation[i][j][k]=0.0;
				}
		}
		for(j=0;j<NDflip;j++){
			for(k=0;k<4;k++){
				Dflip_N_N1[i][j][k]=0;
				Dflip_N_N2[i][j][k]=0;
				Dflip_N_P1[i][j][k]=0;
				Dflip_N_P2[i][j][k]=0;
				Dflip_NAND_Vout[i][j][k]=0.0;
				Dflip_NAND_dissipation[i][j][k]=0.0;
			}
			Dflip_NAND_Vout[i][j][2]=0.0002;
			Dflip_NAND_Vout[i][j][3]=0.0002;
			Dflip_N_N[i][j]=0;
			Dflip_N_P[i][j]=0;
			Dflip_NOT_Vout[i][j]=0.0002;
			Dflip_NOT_dissipation[i][j]=0.0;
			Dflip_dissipation[i][j]=0.0;
		}
	}

	//tree and output_sequence_rearranged data structure for lookup table
	int Ndata = 12;
	vector< vector <vector<int> > > output_tree_ref;
	output_tree_ref.resize(n);
	vector< vector<int> > output_sequence_ref;
	output_sequence_ref.resize(n);
	for(i=0;i<n;i++){
		output_tree_ref.at(i).resize(Ndata/4+2);
		output_tree_ref.at(i).at(0).resize(Ndata);
		output_tree_ref.at(i).at(1).resize(Ndata/2);
		output_tree_ref.at(i).at(2).resize(Ndata/4);
		for(j=3;j<Ndata/4+2;j++){
			output_tree_ref.at(i).at(j).resize(1);
		}
	}

	//generate random sequence
	for(i=0;i<n;i++){
		for(j=0;j<Ndata;j++){
			output_tree_ref[i][0][j]=gsl_rng_uniform_int(rng,2);
		}
	}
	vector< vector<int> > data;
	data.resize(n);
	for(i=0;i<n;i++){
		data[i]=output_tree_ref[i][0];
	}

	if(rank==0){
		printf("input data are:");
		for(j=0;j<data[0].size();j++){
			printf("%d	",data[0][j]);
		}
		printf("\n");
	}

	//fill in look up table for reference
	for(i=0;i<n;i++){
		for(j=0;j<Ndata/2;j++){
			output_tree_ref[i][1][j]=XOR_gate(output_tree_ref[i][0][2*j],output_tree_ref[i][0][2*j+1]);
			output_sequence_ref[i].push_back(output_tree_ref[i][1][j]);
		}
		for(j=0;j<Ndata/4;j++){
			output_tree_ref[i][2][j]=XOR_gate(output_tree_ref[i][1][2*j],output_tree_ref[i][1][2*j+1]);
			output_sequence_ref[i].push_back(output_tree_ref[i][2][j]);
		}
		output_tree_ref[i][3][0]=XOR_gate(output_tree_ref[i][2][0],output_tree_ref[i][2][1]);
		output_sequence_ref[i].push_back(output_tree_ref[i][3][0]);
		for(j=4;j<Ndata/4+2;j++){
			output_tree_ref[i][j][0]=XOR_gate(output_tree_ref[i][j-1][0],output_tree_ref[i][2][j-2]);
			output_sequence_ref[i].push_back(output_tree_ref[i][j][0]);
		}
	}

	if(rank==0){
		printf("output sequence reference is:");
		for(j=0;j<output_sequence_ref[0].size();j++){
			printf("%d	",output_sequence_ref[0][j]);
		}
		printf("\n");
	}

	vector< vector<double> > output_sequence_V_temp, output_sequence_V;
	vector< vector<int> > output_sequence;
	output_sequence_V_temp.resize(n);
	output_sequence_V.resize(n);
	output_sequence.resize(n);

	// VD for Dflip start from 0 and holds the value until changed explicitly
	for(i=0;i<n;i++){
		for(j=0;j<NDflip;j++){
			Dflip_VD[i][j]=0.0;
		}
	}

	int cycle_count=0;
	int Dflip_count=0;

	for(i=0;i<n;i++){
		Dflip_count=0;
		while(data[i].size()>0){
			// send in new data if there is space left on the register
			if(Dflip_count<NDflip){
				XOR_VA[i][0]=undigitalize(data[i][0],Vd);
				XOR_VB[i][0]=undigitalize(data[i][1],Vd);
				XOR_VA[i][1]=undigitalize(data[i][2],Vd);
				XOR_VB[i][1]=undigitalize(data[i][3],Vd);
				data[i].erase(data[i].begin(),data[i].begin()+4);
			}
			// otherwise send data back to ALU and clear up the register
			else{
				XOR_VA[i][0]=output_Dflip_local[i][0];
				XOR_VB[i][0]=output_Dflip_local[i][1];
				XOR_VA[i][1]=output_Dflip_local[i][2];
				XOR_VB[i][1]=output_Dflip_local[i][3];
				Dflip_count=0;
			}
			for(s=0;s<Ntint;s++){
				// first half cycle WE=0, second half cycle WE=1
				for(j=0;j<NDflip;j++){
						if(s<Ntint/2)	Dflip_WE[i][j]=0.01;
						else	Dflip_WE[i][j]=Vd+0.01;
				}
				for(j=0;j<NXOR;j++){
					output_XOR_local[i][j]=XOR_propagation(XOR_N_N1[i][j], XOR_N_N2[i][j], XOR_N_P1[i][j], XOR_N_P2[i][j], XOR_Vout[i][j],XOR_NAND_dissipation[i][j], XOR_dissipation[i][j], XOR_VA[i][j], XOR_VB[i][j], Vd, tint, rng);
				}
				for(j=0;j<NDflip;j++){
					output_Dflip_local[i][j]=Dflip_propagation(Dflip_N_N1[i][j], Dflip_N_N2[i][j], Dflip_N_P1[i][j], Dflip_N_P2[i][j], Dflip_NAND_Vout[i][j], Dflip_NAND_dissipation[i][j], Dflip_N_N[i][j], Dflip_N_P[i][j], Dflip_NOT_Vout[i][j], Dflip_NOT_dissipation[i][j], Dflip_dissipation[i][j], Dflip_VD[i][j], Dflip_WE[i][j], Vd, tint, rng);
				}
				// write on Dflip 0,1 or 2,3 depending on Dflip_count
				Dflip_VD[i][Dflip_count]=output_XOR_local[i][0];
				Dflip_VD[i][Dflip_count+1]=output_XOR_local[i][1];
			}

			// write output sequence
			Dflip_count+=2;
			if(Dflip_count==NDflip){
				for(j=0;j<NDflip;j++){
				output_sequence_V_temp[i].push_back(output_Dflip_local[i][j]);
				}
			}

			if(rank==0 && i==0)	cycle_count+=1;
		}

		//second to last cycle
		XOR_VA[i][0]=output_Dflip_local[i][0];
		XOR_VB[i][0]=output_Dflip_local[i][1];
		XOR_VA[i][1]=output_Dflip_local[i][2];
		XOR_VB[i][1]=output_Dflip_local[i][3];

		for(s=0;s<Ntint;s++){
			for(j=0;j<NDflip;j++){
					if(s<Ntint/2)	Dflip_WE[i][j]=0.01;
					else	Dflip_WE[i][j]=Vd+0.01;
			}
			for(j=0;j<NXOR;j++){
				output_XOR_local[i][j]=XOR_propagation(XOR_N_N1[i][j], XOR_N_N2[i][j], XOR_N_P1[i][j], XOR_N_P2[i][j], XOR_Vout[i][j],XOR_NAND_dissipation[i][j], XOR_dissipation[i][j],XOR_VA[i][j], XOR_VB[i][j], Vd, tint, rng);
			}
			for(j=0;j<NDflip;j++){
				output_Dflip_local[i][j]=Dflip_propagation(Dflip_N_N1[i][j], Dflip_N_N2[i][j], Dflip_N_P1[i][j], Dflip_N_P2[i][j], Dflip_NAND_Vout[i][j], Dflip_NAND_dissipation[i][j], Dflip_N_N[i][j], Dflip_N_P[i][j], Dflip_NOT_Vout[i][j], Dflip_NOT_dissipation[i][j], Dflip_dissipation[i][j], Dflip_VD[i][j], Dflip_WE[i][j], Vd, tint, rng);
			}
			Dflip_VD[i][0]=output_XOR_local[i][0];
			Dflip_VD[i][1]=output_XOR_local[i][1];
		}

		for(j=0;j<2;j++){
			output_sequence_V_temp[i].push_back(output_Dflip_local[i][j]);
		}

		if(rank==0 && i==0)	cycle_count+=1;

		//last cycle
		XOR_VA[i][0]=output_Dflip_local[i][0];
		XOR_VB[i][0]=output_Dflip_local[i][1];

		for(s=0;s<Ntint;s++){
			for(j=0;j<NDflip;j++){
					if(s<Ntint/2)	Dflip_WE[i][j]=0.01;
					else	Dflip_WE[i][j]=Vd+0.01;
			}
			for(j=0;j<NXOR;j++){
				output_XOR_local[i][j]=XOR_propagation(XOR_N_N1[i][j], XOR_N_N2[i][j], XOR_N_P1[i][j], XOR_N_P2[i][j], XOR_Vout[i][j],XOR_NAND_dissipation[i][j], XOR_dissipation[i][j],XOR_VA[i][j], XOR_VB[i][j], Vd, tint, rng);
			}
			for(j=0;j<NDflip;j++){
				output_Dflip_local[i][j]=Dflip_propagation(Dflip_N_N1[i][j], Dflip_N_N2[i][j], Dflip_N_P1[i][j], Dflip_N_P2[i][j], Dflip_NAND_Vout[i][j], Dflip_NAND_dissipation[i][j], Dflip_N_N[i][j], Dflip_N_P[i][j], Dflip_NOT_Vout[i][j], Dflip_NOT_dissipation[i][j], Dflip_dissipation[i][j], Dflip_VD[i][j], Dflip_WE[i][j], Vd, tint, rng);
			}
			Dflip_VD[i][0]=output_XOR_local[i][0];
		}

		output_sequence_V_temp[i].push_back(output_Dflip_local[i][0]);

		if(rank==0){
			if(i==0)	cycle_count+=1;
    			printf("%lf percent completed.\n",100.0*i/n);
  	}
}

	if(rank==0){
		printf("output sequence V temp is:");
		for(j=0;j<output_sequence_V_temp[0].size();j++){
			printf("%lf	",output_sequence_V_temp[0][j]);
		}
		printf("\n");
	}

	//rearrange output sequence
	for(i=0;i<n;i++){
		output_sequence_V[i].push_back(output_sequence_V_temp[i][0]);
		output_sequence_V[i].push_back(output_sequence_V_temp[i][1]);
		output_sequence_V_temp[i][0]=-100.0;
		output_sequence_V_temp[i][1]=-100.0;
		for(j=0;j<Ndata-4;j++){
			if(j%4==2){
				output_sequence_V[i].push_back(output_sequence_V_temp[i][j]);
				output_sequence_V[i].push_back(output_sequence_V_temp[i][j+1]);
				output_sequence_V_temp[i][j]=-100.0;
				output_sequence_V_temp[i][j+1]=-100.0;
			}
		}
		output_sequence_V[i].push_back(output_sequence_V_temp[i][4]);
		output_sequence_V_temp[i][4]=-100.0;
		for(j=5;j<Ndata-1;j++){
			if(j%4==1){
				output_sequence_V[i].push_back(output_sequence_V_temp[i][j]);
				output_sequence_V_temp[i][j]=-100.0;
			}
		}
		for(j=0;j<output_sequence_V_temp[i].size();j++){
				if(output_sequence_V_temp[i][j]>-99.0){
				output_sequence_V[i].push_back(output_sequence_V_temp[i][j]);
				}
		}
	}

	if(rank==0){
		printf("output sequence V is:");
		for(j=0;j<output_sequence_V[0].size();j++){
			printf("%lf	",output_sequence_V[0][j]);
		}
		printf("\n");
	}

	//covert output to binary and compare to reference
	vector<int> error_final_local(n),error_count_local(n);
	vector<int> error_final_global(Nw),error_count_global(Nw);
	double error_avg=0.0;
	vector<double> dissipation_tot_local(n),dissipation_tot_global(Nw);
	double dissipation_tot_avg=0.0;

	for(i=0;i<n;i++){
		error_count_local[i]=0;
		error_final_local[i]=0;
		dissipation_tot_local[i]=0.0;
		for(j=0;j<Ndata-1;j++){
			output_sequence[i].push_back(digitalize(output_sequence_V[i][j],Vd));
			if(output_sequence[i][j]!=output_sequence_ref[i][j]){
				output_sequence[i][j]=-1;
				error_count_local[i]+=1;
				if(j==Ndata-2)	error_final_local[i]=1;
			}
		}
		for(j=0;j<NXOR;j++){
			dissipation_tot_local[i]+=XOR_dissipation[i][j];
		}
		for(j=0;j<NDflip;j++){
			dissipation_tot_local[i]+=Dflip_dissipation[i][j];
		}
	}
	MPI_Gather(&error_final_local[0],n,MPI_INT,&error_final_global[0],n,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&error_count_local[0],n,MPI_INT,&error_count_global[0],n,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&dissipation_tot_local[0],n,MPI_DOUBLE,&dissipation_tot_global[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if(rank==0){
		printf("output sequence is:");
		for(j=0;j<output_sequence[0].size();j++){
			printf("%d	",output_sequence[0][j]);
		}
		printf("\n");
		printf("error final is:%d\n",error_final_local[0]);
		printf("error count is:%d\n",error_count_local[0]);
		printf("total dissipation is:%lf\n",dissipation_tot_local[0]);

		for(i=0;i<Nw;i++){
			error_avg+=1.0/Nw*error_final_global[i];
			dissipation_tot_avg+=1.0/Nw*dissipation_tot_global[i];
		}
		printf("overall average error rate is %lf, total dissipation is %lf\n",error_avg,dissipation_tot_avg);
	}

	//printout output_sequence
	vector< vector<int> > output_local(Ndata-1, vector<int>(n));
	vector< vector<int> > output_global(Ndata-1, vector<int>(Nw));
	vector< vector<double> > output_V_local(Ndata-1, vector<double>(n));
	vector< vector<double> > output_V_global(Ndata-1, vector<double>(Nw));
	for(j=0;j<Ndata-1;j++){
		for(i=0;i<n;i++){
			output_local[j][i]=output_sequence[i][j];
			output_V_local[j][i]=output_sequence_V[i][j];
		}
	}
	for(j=0;j<Ndata-1;j++){
		MPI_Gather(&output_local[j][0],n,MPI_INT,&output_global[j][0],n,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(&output_V_local[j][0],n,MPI_DOUBLE,&output_V_global[j][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	}

	if(rank==0){
		result=fopen("/global/scratch/chloegao/traj_device_Vd=8_2XOR_12D_T=1e7","w");
		for(i=0;i<Nw;i++){
			fprintf(result,"%d	%d	%lf	",error_final_global[i],error_count_global[i],dissipation_tot_global[i]);
			for(j=0;j<Ndata-1;j++){
				fprintf(result,"%d	",output_global[j][i]);
			}
			for(j=0;j<Ndata-1;j++){
				fprintf(result,"%lf	",output_V_global[j][i]);
			}
			fprintf(result,"\n");
		}
		fclose(result);
	}

	MPI_Finalize();
	return 0;
}
