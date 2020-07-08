#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <cmath>
#include "mpi.h"
using namespace std;
	  
	  double dt = 0.0001; 
	  double sd = 0.01; 
	  int T = 1;
	  const int N = 10000; 
	  int HN = 5000;  
	  const int M = 300000;  ///          on each proc thus 30 means 300 000
	  const int Nc= 30;     
	  const int p = 30;    

void BM_SK(const double, double [N], double [N+1]);

void AsmG_P(double [N+1],double [N+1],double [N], double [N+1], double [M], double [N+1], double [2], double [2], const double);

void Naive_Eu(double W[N+1],double dW[N], double Y[N+1],  double a[2], double b[2]);
void Eu(double [N+1],double [N], double [N+1], double [2], double [2]);

void P_Eu_Fine(int, double [Nc+1], double [Nc+1], double [N+1],double [N], double [N+1], double [2], double [2], int);
void P_Eu_Coarse(double [N+1], double [Nc+1], double, double [2], double [2], int);
void P_er_pro(double [Nc+1], double [N+1], double [Nc+1], double [Nc+1], int, double, double [2], double [2], double [N+1]);

void P_Naive_Eu_Coarse(double [N+1], double [Nc+1], double, double [2], double [2], int);
void P_Naive_Eu_Fine(double [Nc],int, double [Nc+1], double [Nc+1], double [N+1],double [N], double [N+1], double [2], double [2], int);
void P_Naive_Eu_Pro_MinPaths(double [Nc+1],double [Nc]);
void P_Naive_Eu_er_pro(double [Nc+1], double [N+1], double [Nc+1], double [Nc+1], int, double, double [2], double [2], double [N+1]);

void Ex_sol(double [N+1],double [N], double [N+1], double [2], double [2]);

      
int main(int argc, char **argv)
{  

	  
	//-----------------------------------------------------------------------------------
	//************************ Initialis      *******************************************
	//-----------------------------------------------------------------------------------
    // 1- ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  int pp;    
	  int myrank;
	  int source;
	  int dest = 0;
	  int tag = 0;
	  MPI_Status status;
	  MPI_Request myRequest;

	  MPI_Init(&argc, &argv);
	  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	  MPI_Comm_size(MPI_COMM_WORLD, &pp);
	
	// 2- Stc Pr   +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++		
	  double X[N+1];   
	  double Y[N+1];   	  
	  double W[N+1];  
	  double dW[N];  
      double YT[M]; 
        
	// 3- param     +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
	  double a[2] = {4, -2};  
	  double b[2] = {2, 1}; 
	  //double b_2[2] = {0,0}; 
	  //a[0]=1; 
	  //a[1]= -2;
	  
   	// 4- NRG  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  
	  srand((unsigned)time(0));
      const double scale = 1/double(RAND_MAX);

  	// 5- Env ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  double lamda[2] = {1, 1};           
	  double dtz[N];
	  double tz[N+1];
	  int size_tz = 0;
	  tz[0] = 0;

	// 6- Par   ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 	  int NNc = 10;           
	  double Dtc = 0.025;                             
	  double Sc[Nc+1];                      // size Nc
	  double Yc[Nc+1];                      // size Nc
	  double Yc2[Nc+1];
	  
	         	
	// 7- P-Naive-Eu +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	  double Min[Nc];


	//-----------------------------------------------------------------------------------
	//************************  Main Program  *******************************************
	//-----------------------------------------------------------------------------------
     
	  //for (int ik =0; ik < NN; ik++){
		  //b_2[0] = b_2[0] + 0.2;
		  //b_2[1] = b_2[1] + 0.2;
		  //b[0] = sqrt(b_2[0]); 
		  //b[1] = sqrt(b_2[1]);
	
	  //for (int jj=0; jj< M; jj++){  

	  Y[0] = 1;
	  Yc[0] = Y[0];
	  Yc2[0] = Y[0];
	  X[0] = Y[0];
	  W[0] = 0;

	  for (int i = 0; i<Nc+1; i++){
		  Sc[i] =0;                         //used in the P_N_E
	  }

	   BM_SK(scale, dW, W);

	  //Y[0] = 0; 
	   //   rr = 0;
		//  while (rr ==0)
		//  {
		//	  rr =  rand()*scale;
		//  }
		//   Y[0] = -0.5*log(rr); 


	  //Naive_Eu(W, dW, Y, a, b);
	  //Ex_sol(W, dW, Y, a, b);
	  Eu(W, dW, Y, a, b);

	  if (myrank == 0 ){
		  string fname;
		  string filename;			  
		  char str[33];
		  sprintf (str, "%d", M);
		  filename = str;
		  fname = "serial" + filename ;
		  ofstream my1file;
		  my1file.open (fname.c_str());

	  for (int i=0; i<N+1; i++) {  
	       my1file << Y[i] << endl;
	  }
	  my1file.close();
	  }


	  P_Eu_Coarse(W, Yc, Dtc, a, b, NNc);	  
	  // P_Naive_Eu_Coarse(W, Yc, Dtc, a, b, NNc);

	   if (myrank == 0){
			int i0c = myrank;
		    // for (ic=0; ic < Nc; ic++){			  
			P_Eu_Fine(i0c, Yc, Yc2, W, dW, Y, a, b, NNc);
	   }
	  
	  for (int iu = 0; iu<4; iu++){
 		  if (myrank ==0){
			  for (source=1; source<pp; source++){
				 int ipp = source;   // ip = Nc
				  MPI_Recv(&Yc2[ipp+1], 1, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);	
			  }
		  }
		  if (myrank != 0){
			  int icc = myrank;
			  // for (ic=0; ic < Nc; ic++){			  
			  P_Eu_Fine(icc, Yc, Yc2, W, dW, Y, a, b, NNc);
			  //P_Naive_Eu_Fine(Min, ic, Yc, Yc2, W, dW, Y, a, b, NNc);
			  MPI_Send(&Yc2[icc+1], 1, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			  }
		 	  
		  
		  if (myrank ==0){
		  //P_Naive_Eu_Pro_MinPaths(Yc2, Min);
		  P_er_pro(Yc, Y, Yc2, Sc, NNc, Dtc, a, b, W);
		   }
		   
		   MPI_Bcast(Yc, (Nc+1), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		   MPI_Barrier(MPI_COMM_WORLD);
	  }
     
  
      
	  if (myrank  != 0){
		  int imc = myrank;
		  // for (ic=0; ic < Nc; ic++){	
		  P_Eu_Fine(imc, Yc, Yc2, W, dW, Y, a, b, NNc);
		  //P_Naive_Eu_Fine(Min, ic, Yc, Yc2, W, dW, Y, a, b, NNc);
	  }
	  
	   if (myrank != 0){		   	  	 
			  //................... OUTPUT ......................//
			  int ipi = myrank;
			  string fname;
			  string filename;			  
				  char str[33];
				  sprintf (str, "%d", myrank);
				  filename = str;
			 	  fname = filename;
				  ofstream myfile;
				  myfile.open (fname.c_str());
				  				   			  
				  for(int i=((ipi*NNc)+1); i<((ipi+1)*NNc+1); i++){
					  myfile  <<"proc  "<< ipi << "  Y=  "<< Y[i] << endl;
                    }

			  myfile.close();
			  }

	    if (myrank ==0){
			  for (source=1; source<pp; source++){
				   int ii0 = source;   // ii = Nc
				   MPI_Recv(&Y[(ii0*NNc)+1], NNc, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);	
			          }
		            }
		  if (myrank != 0){
			  int ii0 = myrank;
			  MPI_Send(&Y[(ii0*NNc)+1], NNc, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			  }

	 
		  
          if (myrank ==0){	
		  string fname;
		  string filename;			  
		  char str[33];
		  sprintf (str, "%d", M);
		  filename = str;
		  fname = "pl" + filename ;
		  ofstream my2file;
		  my2file.open (fname.c_str());	
		  for (int i=0; i<N+1; i++) {            
		    my2file << Y[i] << endl;
		  }
		  my2file.close();
	   }

  // }
   

	  //cout << myrank << endl; 
	  
			   //................................................//				
	  //YT[jj] = Y[N-1];
	   // for (int ii=0; ii<N; ii++){
		//  myfile <<  tz[ii] << "     " << Y[ii] << "    " << X[ii] << endl;	   // }    
	  //EYT[ik] = 0;
      //for (int i=0; i<M; i++){
		  //EYT[ik] = EYT[ik] + YT[i];	  //}
	  //EYT[ik] = EYT[ik]/M;	  
	  //for (int i=0; i<NN; i++){
		//  myfile <<  EYT[i] << endl;	  //}
	  //char f;
	  //cin >> f;
	  //cout << myrank << "done" << endl;
      MPI_Finalize();
  	  return 0;      
}

	void BM_SK(const double scale, double dW[N], double W[N+1]){
	  double U1; 
	  double U2;
	  double V1;
	  double V2;
	  double S;
	  double X1;
	  double X2;

      for(int i=0; i< HN; i++){
		  S = 2.0;
		  while (S >= 1 || S ==0)
		  {
			  U1 =  rand()*scale;
			  U2 =  rand()*scale;
			   
			  V1 = 2*U1-1;
			  V2 = 2*U2-1;
			  S = (V1*V1)+(V2*V2);
			   
				  if (S < 1 && S >0)
				  {
					  X1 = V1 * sqrt ((-2*log(S))/S);  
					  X2 = V2 * sqrt ((-2*log(S))/S);					  
				  }
		  }
		  
		  dW[(2*i)] = sd *X1;
		  dW[(2*i)+1] = sd *X2; 
	  }
  
	  for (int i=0; i<N; i++){
	  W[i+1] = W[i] + dW[i];
	  }
	   
	}



	  void AsmG_P(double tz[N+1], double W[N+1],double dW[N], double Y[N+1], double YT[M], double X[N+1], double a[2], double b[2], const double scale){
      double R;        // parameter useless inside this method
	  double Un;	   // parameter useless inside this method
      int i=0;
	  double t_iplus1;
	  t_iplus1 = dt;

	  		int j=1;
			int k = 1;
			int jz = 1;
	
	  while ( i < N){
		  j = j-k;
		  
			while (t_iplus1 < tz[jz]){

				dW[i] = (a[j] *dt) + (b[j] * dW[i]); // this is Delta_W step (ii-1) of the approximation of Y[i+1]
				
				Un = 0;
				while (Un == 0 || Un == 1){
					Un = rand()*scale;
				}
				R = (-dW[i] /2)+( (sqrt((dW[i]*dW[i]) - (2*dt* log(Un))))/2);
				
				//mmmyfile <<  Y[0] << "     " << R <<  endl;
				R = R - Y[i];
				if( R < 0){
					R = 0;
				}
				Y[i+1] = Y[i] +  dW[i] + R; // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
				X[i+1] = X[0] + (a[j] *i*dt) + (b[j]*W[i+1]); //i*dt = t 

				i = i+1;
				t_iplus1 = i*dt;
			}

				k = k*(-1);
				jz = jz + 1;				
			}			
	  }
	
	  
	  

	  void Naive_Eu(double W[N+1],double dW[N], double Y[N+1],  double a[2], double b[2]){
      double R;
	  double ddw;
	  int j=0;
	  for (int i=0; i < N; i++){
				ddw = (a[j] *dt) + (b[j]*(Y[i]) * dW[i]); // this is Delta_W step (ii-1) of the approximation of Y[i+1]							
				//mmmyfile <<  Y[0] << "     " << R <<  endl;
				R = -Y[i] -ddw;
				if( R < 0){
					R = 0;
				}
				Y[i+1] = Y[i] +ddw + R; // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
			}			
	  }

	  void P_Naive_Eu_Coarse(double W[N+1], double Yc[Nc+1], double Dtc, double a[2], double b[2], int NNc){
	      double R;
		  double ddw;
		  int j= 0;
		  for(int i=0; i< Nc; i++){			  	
				ddw = (a[j]*Dtc) + (b[j]*(Yc[i])* (W[(i+1)*NNc]-W[i*NNc])); // this is Delta_W step (ii-1) of the approximation of Y[i+1]				
				R = -Yc[i] -ddw;
				if( R < 0){
					R = 0;
				}
				Yc[i+1] = Yc[i] +  ddw + R;				
		  }
	  }

	  void P_Naive_Eu_Fine(double Min[Nc], int ic, double Yc[Nc+1], double Yc2[Nc+1], double W[N+1],double dW[N], double Y[N+1],  double a[2], double b[2], int NNc){
	  double R;
	  int j=0;
	  double ddw;
	  Y[ic*NNc] = Yc[ic];
	  Min[ic] = 0;
		   for(int i=(ic*NNc); i< (ic+1)*NNc; i++){	
			   ddw = (a[j] *dt) + (b[j]*(Y[i]) *dW[i]); // this is Delta_W step (ii-1) of the approximation of Y[i+1]				
						
				//mmmyfile <<  Y[0] << "     " << R <<  endl;
				R = -Y[i] -ddw;
				if( R < 0){
					R = 0;
				}
				Y[i+1] = Y[i] +  ddw + R; // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
				Min[ic] = Min[ic] + R;
			}		
		   		   Yc2[ic+1] = Y[(ic+1)*NNc];
	  }

	  void P_Naive_Eu_Pro_MinPaths(double Yc2[Nc+1],double Min[Nc]){
		  double R1;
		  for(int i=1; i<Nc; i++){
			  	R1 = Min[i-1] - Min[i];
				if( R1 < 0){
					R1 = 0;
				}
				Min[i] = Min[i] + R1;
				////////////Yc2[i+1] = Yc2[i+1] ____??????????????????///// + R1; 
		  }
	  }

	  	 void P_Naive_Eu_er_pro(double Yc[Nc+1], double Y[N+1], double Yc2[Nc+1], double Sc[Nc+1], int NNc, double Dtc, double a[2], double b[2], double W[N+1]){
		 int j= 0;
		 double R;
		 double ddw;
		 for (int i=1; i<Nc; i++){
			 Sc[i] =  Sc[i] +   (Yc2[i] - Yc[i]);
		 }
		  for(int i=0; i< (Nc-1); i++){			  	
				ddw = (a[j]*Dtc) + (b[j]*(Yc[i])*(W[(i+1)*NNc]-W[i*NNc]))+ Sc[i+1]; // this is Delta_W step (ii-1) of the approximation of Y[i+1]				
				R = -Yc[i] -ddw;
				if( R < 0){
					R = 0;
				}
				Yc[i+1] = Yc[i] +  ddw + R;				
		  }
		 }



	  void Eu(double W[N+1],double dW[N], double Y[N+1], double a[2], double b[2]){ 
		   int j= 0;
	       for(int i=0; i< N; i++){		       	 
		 						
		 Y[i+1] = Y[i] +  (a[j]*dt) + (b[j]* dW[i]); // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
				// X[i+1] = X[0] + (a[j] *i*dt) + (b[j]*W[i+1]); //i*dt = t 
				}			
	  }
	  void Ex_sol(double W[N+1],double dW[N], double Y[N+1], double a[2], double b[2]){
		   int j= 0;
	       for(int i=0; i< N; i++){		       	 
		 		Y[i+1] = Y[i]*exp(((a[j]-(b[j]*b[j]/2))*dt) + (b[j] * dW[i]));  // exact solution
				}			
	  }

      void P_Eu_Fine(int ic, double Yc[Nc+1], double Yc2[Nc+1], double W[N+1],double dW[N], double Y[N+1], double a[2], double b[2], int NNc){
		   int j= 0;
		   Y[ic*NNc] = Yc[ic];
		   for(int i=(ic*NNc); i< (ic+1)*NNc; i++){		       	 
		 						
				Y[i+1] = Y[i] +  (a[j]*dt) + (b[j]* dW[i]); // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
				// X[i+1] = X[0] + (a[j] *i*dt) + (b[j]*W[i+1]); //i*dt = t 
				}	
		   Yc2[ic+1] = Y[(ic+1)*NNc];
	  }

	 void P_Eu_Coarse(double W[N+1], double Yc[Nc+1], double Dtc, double a[2], double b[2], int NNc){
	      int j= 0;
		  for(int i=0; i< Nc; i++){
			Yc[i+1] = Yc[i] +  (a[j]*Dtc) + (b[j]* (W[(i+1)*NNc]-W[i*NNc])); // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
			// X[i+1] = X[0] + (a[j] *i*dt) + (b[j]*W[i+1]); //i*dt = t 
		  }
	  }

	 void P_er_pro(double Yc[Nc+1], double Y[N+1], double Yc2[Nc+1], double Sc[Nc+1], int NNc, double Dtc, double a[2], double b[2], double W[N+1]){
		 for (int i=1; i< (Nc+1); i++){
			 Sc[i] = Sc[i] + (Yc2[i] - Yc[i]);
		 }
		 int j= 0;
		  for(int i=0; i<Nc; i++){
				Yc[i+1] = Yc[i] +  (a[j]*Dtc) + (b[j]*(W[(i+1)*NNc]-W[i*NNc]))+ Sc[i+1]; // dW[i] = W[i+1] - W[i]; this is Delta_W step (ii-2) of the approximation of Y[i+1]
				// X[i+1] = X[0] + (a[j] *i*dt) + (b[j]*W[i+1]); //i*dt = t 
		  }


	 }


