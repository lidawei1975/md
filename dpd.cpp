
//#define SSL
//#define MYMPI	

#define SMALL 0.001f
#define TYPE            8	//total types of atoms +1
#define TYPEPOLY        2	//total types of protein +1

//#define RELAX
//#define PRESSURE
//#define PRESSURECOUP
//#define ANISOTROPIC
#define TEMPERATURECOUP
//#define MYDEBUG

//#define LANGEVIN
#define HARDCORE
#define HYDROPHOBIC
#define INTERACTION14
#define HBOND
#define BOND
#define BONDANGLE
#define TORSION
#define OUTOFPLANE
//#define SHAKE


#include <cmath>
#include <fstream>
#include <cstring>


#ifdef MYMPI
#include "mpi.h"
#endif				 

#ifdef SSL
#include "xmmintrin.h"
#endif

using namespace std;

clock_t t1,t2;

float *rx,*ry,*rz,*rvx,*rvy,*rvz,*rax,*ray,*raz,**rold,*f;
int   ** cell;

int monomernumber;
int minimizer,Min_Step;
float region[4],regionH[4];
float volume,svolume,ssvolume;
float rCut,rCut_second,rrCut,rrCut_second,deltaT,sqrtdeltaT,ddeltaTH;
float Langevin,beta,lamad,density,pi,temperature,timeNow,coup_temperature,coup_pressure;
float vMag,vSum,vvSum;

#ifdef PRESSURECOUP
#ifdef ANISOTROPIC
float miu[4];
#else
float miu;
#endif
#endif

float fa,sinfa0,cosfa0;
float kinEnergy,potEnergy,bondEnergy;  //these for calculate the energy, velocity
float hEnergy,hydroEnergy,angleEnergy,torsionEnergy,hNumber,outofplaneNumber;
float vVal,uVal,hVal,hydroVal,potVal,bondVal,angleVal,torsionVal,hnumberVal,outofplaneVal;
float sKinEnergy,ssKinEnergy,pressure,spressure,sspressure;
long  seed=-5;
float virSum[4][4];
float kinSum[4][4];
float svirSum[4][4];
float skinSum[4][4];

int   BuildListReturn;
int   LeapfrogStepReturn;
int   ComputeBondForceReturn;
int   beginnumber;
int   dim;
int   version;
int   moreCycles,nAtom;
int   stepAvg,stepCount,stepEquil,stepLimit;
int   *type,*type2,*inPoly;
char  use[500];
int   byteread;  //byte_read in structure file
ofstream fout;  //main file for output informations.
FILE  *fp;       //used for read from structures file.


int   *cellList,cells[4],cells_second[4];                  //this line for cell method
float dispHi,rNebrShell,vvMax;                 //neighbor list method
float dispHi_second,rNebrShell_second;                 //neighbor list method
int   *nebrTab1,*nebrTab2,*nebrTab3,*nebrTab4;                      //neighbor list method
int   nebrNow,nebrTabFac,nebrTabLen,nebrTabLen2,nebrTabMax; //neighbor list method
int   *nebrTab1_second,*nebrTab2_second,*nebrTab3_second,*nebrTab4_second;                      //neighbor list method
int   nebrNow_second,nebrTabLen_second,nebrTabLen2_second; //neighbor list method 


//read the force field file and topological file to fill these values !!
int   TypePoly;
float constraintPrec;
int   nChain[TYPEPOLY],ChainLength[TYPEPOLY],ChainBegin[TYPEPOLY];
int   AtomNumber[TYPEPOLY];
int   BondNumber[TYPEPOLY];
int   BondAngleNumber[TYPEPOLY];
int   BondTorsionNumber[TYPEPOLY];
int   Interaction14Number[TYPEPOLY];
int   OutofPlaneNumber[TYPEPOLY];
int   ShakeNumber[TYPEPOLY];
int   ForceNumber;
int   TypeAtom;

float mm[TYPE];

int   *BondIndex1[TYPEPOLY];
int   *BondIndex2[TYPEPOLY];
float *BondLength[TYPEPOLY];
float *BondStrenth[TYPEPOLY];

int   *BondAngleIndex1[TYPEPOLY];
int   *BondAngleIndex2[TYPEPOLY];
int   *BondAngleIndex3[TYPEPOLY];
float *BondAngle[TYPEPOLY];
float *BondAngleStrenth[TYPEPOLY];

int   *BondTorsionIndex1[TYPEPOLY];
int   *BondTorsionIndex2[TYPEPOLY];
int   *BondTorsionIndex3[TYPEPOLY];
int   *BondTorsionIndex4[TYPEPOLY];
int   *BondTorsionType[TYPEPOLY];
float *BondTorsionStrenth[TYPEPOLY];

int   *OutofPlaneIndex1[TYPEPOLY];
int   *OutofPlaneIndex2[TYPEPOLY];
int   *OutofPlaneIndex3[TYPEPOLY];
int   *OutofPlaneIndex4[TYPEPOLY];
float *OutofPlaneStrenth[TYPEPOLY];

int   *Interaction14Index1[TYPEPOLY];
int   *Interaction14Index2[TYPEPOLY];

int   *ShakeIndex1[TYPEPOLY];
int   *ShakeIndex2[TYPEPOLY];
float *ShakeLength[TYPEPOLY];

float vmForceMatrix[TYPE][TYPE];
float elForceMatrix[TYPE][TYPE];
float vmDeltaMatrix[TYPE][TYPE];
float vmForceMatrix14[TYPE][TYPE];
float vmDeltaMatrix14[TYPE][TYPE];
float H_Force,Hydro_Force;

int   **ExcludeMatrix[TYPEPOLY];
//////////////////////////////////////////////////////////////////////////

//For script
char  initfilename[300];
char  messfilename[300];
char  maindirectory[300];

#ifdef SSL
int   rr_to_fcVal(float *, float *, float *, float *);
#endif

int   DPDmain();
int   ProcessError(int);
void  AllocateMemory();
void  ReleaseMemory();
void  AccumProps(int);
void  ApplyBoundaryCond();
int   Leapfrog();
int   shake();
int   ComputeLangevin();
int   ComputeForce_hardcore();
int   ComputeForce_hydrophobic();
int   ComputeWallForce();
int   ComputeBondForce();
int   BuildNebrList();
int   BuildNebrList_second();
void  EvalProps();
void  EvalBuildNeighbor();
void  EvalBuildNeighbor_second();
void  OutPut(int Index);
void  swap(int p1, int p2);
void  AddPolymer();
void  AdjustZ();
float Integrate(float *f,int nf);  //math
double ran(long *idum);
void  Gettension(int flag, float f[4], float dr[4], int j1, int j2);
void  Outputdata();
void  Receivecommand();
int   Setparameter();
void  myprintf(char *ps);

#ifdef MYMPI
int nprocs,cube;
int ncube,nlocal;
int   nchainstart[10],nchainend[10];
float *raxx,*rayy,*razz;
void FoldForce(void);
void Expandr(void);
#endif

int node,nstart,nend;


main(int argc, char * argv[])
{

#ifdef MYMPI
	MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&node);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
      
	ncube = 0;
	if(nprocs==2)
		ncube=1;
	else
	{
		sprintf(use,"error! \n");
		myprintf(use);
		MPI_Finalize();
		return 0;
	}
#else
	node=0;	
#endif

    if (node == 0)
        printf(use, "My DPD program, version 2.51 begin \r\n");

    if (node == 0)
        printf("Perfomr mission in current directory\r\n");
    sprintf(initfilename, "./dpd.ini");
    sprintf(messfilename, "./mess.txt");
    sprintf(maindirectory, "./");
    DPDmain();

#ifdef MYMPI
	MPI_Finalize();
#endif
	sprintf(use,"\r\nAll missions are finished\r\n");
	myprintf(use);
	return 0;
}


int DPDmain()
{
	if(Setparameter()==0) 
	{
#ifdef MYMPI
		MPI_Finalize();
#endif
		return 0;	
	}

	float h,uOld,fmax;
	int i,n,k;
	int times,timeh,timem;
	if(minimizer==1)
	{
		
		rold=new float* [nAtom+1];
		for(k=1;k<=3;k++)
		{
			rold[k]=new float[nAtom+1];
		}
		f=new float[nAtom+1];

		h=0.01f;
		uOld=1000000000.0f;
		AccumProps(0);
		printf("Now begin minimize, total %d steps\r\n",Min_Step);
		for(i=0;i<=Min_Step;i++)
		{
			BuildNebrList();
			ComputeBondForce();
			EvalProps();
			AccumProps(1);

			if(uVal+bondVal<uOld)
			{
				h=1.2f*h;
				uOld=uVal+bondVal;
				if(i%10==0 && i!=0) printf("step is %d, f is %f, potential is %f\r\n",i,fmax,uOld);
			}
			else
			{
				h=0.2f*h;
				for(n=1;n<=nAtom;n++)
				{
					rx[n]=rold[1][n];
					ry[n]=rold[2][n];
					rz[n]=rold[3][n];
				}
				BuildNebrList();
				ComputeBondForce();
				EvalProps();
				AccumProps(1);
				uOld=uVal+bondVal;
			}

			fmax=0.0f;
			for(n=1;n<=nAtom;n++)
			{
				if(n==82)
					n=82;
				f[n]=(float)sqrt(rax[n]*rax[n]+ray[n]*ray[n]+raz[n]*raz[n]);
				if(f[n]>fmax)
					fmax=f[n];
			}

			for(n=1;n<=nAtom;n++)
			{
				rold[1][n]=rx[n];
				rold[2][n]=ry[n];
				rold[3][n]=rz[n];
				rx[n]+=h*rax[n]/fmax;
				ry[n]+=h*ray[n]/fmax;
				rz[n]+=h*raz[n]/fmax;
			}

			ApplyBoundaryCond();
		}
		OutPut(-1);
		for(k=1;k<=3;k++)
		{
			delete rold[k];
		}
		delete rold;
		delete f;	
	}


	sprintf(use,"Doing dynamics!\r\n");
	myprintf(use);

	sprintf(use,"time  Kin     Pot    h     hNumber  hydro    bond    angle  torsion outplane");
	myprintf(use);
	
	if(node==0) printf("\r\n");
	sprintf(use,"   vir1    kin1    vir2    kin2   vir3     kin3");
	
	fout<<use;
	fout<<"\r\n";
	fout.flush();

	AccumProps(0);
	stepCount=0;
	nebrNow=1;
	nebrNow_second=1;
	t1=clock();
	while(moreCycles)
	{
		stepCount=stepCount+1;
		timeNow=stepCount*deltaT;

		if(nebrNow)
		{
			BuildListReturn=BuildNebrList();
			if(ProcessError(2)==0) return 0;
		}

		if(nebrNow_second)
		{
			BuildListReturn=BuildNebrList_second();
			if(ProcessError(2)==0) return 0;
		}

		for(n=1;n<=nAtom;n++)
		{
		rax[n]=0.0f;
		ray[n]=0.0f;
		raz[n]=0.0f;
		}

		uVal=0.0f;
		hVal=0.0f;
		hnumberVal=0.0f;
		hydroVal=0.0f;
		bondVal=0.0f;
		angleVal=0.0f;
		torsionVal=0.0f;
		outofplaneVal=0.0f;
	
	    ComputeForce_hardcore();	  //hardcore and H-bond potential
				

#ifdef HYDROPHOBIC
        ComputeForce_hydrophobic();
#endif

		ComputeBondForceReturn=ComputeBondForce();

#ifdef LANGEVIN
		ComputeLangevin();
#endif


#ifdef MYMPI
		FoldForce();
#endif

		if(ProcessError(3)==0) return 0;

		AccumProps(1);

		LeapfrogStepReturn=Leapfrog(); 
		if(ProcessError(4)==0) return 0;

#ifdef SHAKE
		shake();
#endif
		ApplyBoundaryCond();

#ifdef MYMPI
		Expandr();
#endif

		EvalBuildNeighbor();
		EvalBuildNeighbor_second();

		if(stepCount%stepAvg == 0)
		{	
			AccumProps(2);
			Outputdata();
			if(stepLimit!=0 && stepCount>=stepLimit)
				moreCycles=0;
			AccumProps(0);
			Receivecommand();
		}
	}
	sprintf(use,"cycle finished, release memory \r\n");
	myprintf(use);
	ReleaseMemory();
	sprintf(use,"successfully finished!");
	myprintf(use);
	t2=clock();
	times=(t2-t1)/CLOCKS_PER_SEC;
	timeh=times/3600;
	timem=(times-3600*timeh)/60;
	times=times-3600*timeh-60*timem;
	sprintf(use," **** The total time is %dhours %dmimutes and %dseconds!\r\n\r\n\r\n",timeh,timem,times);
	myprintf(use);
	fout.close();
	return 1;
}

int Leapfrog()
{
	int n;
	for(n=nstart;n<=nend;n++)
	{

		rax[n]*=mm[type[n]];
		ray[n]*=mm[type[n]];
		raz[n]*=mm[type[n]];

#ifdef RELAX 
		rvx[n]=rvy[n]=rvz[n]=0.0f;
#endif
		rvx[n]+=deltaT*rax[n];
		rvy[n]+=deltaT*ray[n];
		rvz[n]+=deltaT*raz[n];

#ifdef TEMPERATURECOUP
		rvx[n]*=lamad;
		rvy[n]*=lamad;
		rvz[n]*=lamad;
#endif
		rx[n] +=deltaT*rvx[n];
		ry[n] +=deltaT*rvy[n];
		rz[n] +=deltaT*rvz[n];

	}


#ifdef PRESSURECOUP	
#ifdef ANISOTROPIC
	for(int k=1;k<=3;k++)
	{
		region[k]*=miu[k];
		regionH[k]=region[k]/2;
	}
	step_tension*=miu[3];
#else
	for(int k=1;k<=3;k++)
	{
		region[k]*=miu;
		regionH[k]=region[k]/2;
	}
	step_tension*=miu;
#endif
	maxtensionstep=1+(int)(region[3]/step_tension);
#endif

	return 0;
}

int  shake()
{
	float dr[4];
	float cDev,cDevR,cDevV,g,ga;
	int   changed,i,m,maxCycle,j1,j2,n,nn;
	int   cycleR,cycleV;

	maxCycle=200;
	cDevR=cDevV=0.0f;

	for(nn=1;nn<=TypePoly;nn++)
	{  
		for(n=1;n<=nChain[nn];n++)
		{	 
			m=ChainBegin[nn]+(n-1)*ChainLength[nn];
			cycleR=0;
			changed=1;
			while(cycleR<maxCycle && changed)
			{
				cycleR++;
				changed=0;
				cDev=0.0f;
				for(i=1;i<=ShakeNumber[nn];i++)
				{
					j1=m+ShakeIndex1[nn][i];
					j2=m+ShakeIndex2[nn][i];  

					dr[1]=rx[j1]-rx[j2];
					if(dr[1]>regionH[1])
						dr[1]-=region[1];
					else if(dr[1]<-regionH[1])
						dr[1]+=region[1];
					
					dr[2]=ry[j1]-ry[j2];
					if(dr[2]>regionH[2])
						dr[2]-=region[2];
					else if(dr[2]<-regionH[2])
						dr[2]+=region[2];
					
					dr[3]=rz[j1]-rz[j2];
					if(dr[3]>regionH[3])
						dr[3]-=region[3];
					else if(dr[3]<-regionH[3])
						dr[3]+=region[3];

					g=(dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3]-ShakeLength[nn][i])/(4.0f*ShakeLength[nn][i]);
					ga=(float)fabs(g);
					if(ga>cDev)
						cDev=ga;
					if(ga>constraintPrec)
					{
						changed=1;
						rx[j1]-=g*dr[1];
						rx[j2]+=g*dr[1];
						ry[j1]-=g*dr[2];
						ry[j2]+=g*dr[2];
						rz[j1]-=g*dr[3];
						rz[j2]+=g*dr[3]; 
					}
				}
			}

			if(cycleR>=maxCycle) 
				printf("shake R mistake at step %d!! and cDev is %f \n",stepCount,cDev);

			if(cDev>cDevR)
				cDevR=cDev;
			
			cycleV=0;
			changed=1;

			while(cycleV<maxCycle && changed)
			{
				cycleV++;
				changed=0;
				cDev=0.0f;
				for(i=1;i<=ShakeNumber[nn];i++)
				{
					j1=m+ShakeIndex1[nn][i];
					j2=m+ShakeIndex2[nn][i];

					dr[1]=rx[j1]-rx[j2];
					if(dr[1]>regionH[1])
						dr[1]-=region[1];
					else if(dr[1]<-regionH[1])
						dr[1]+=region[1];
					
					dr[2]=ry[j1]-ry[j2];
					if(dr[2]>regionH[2])
						dr[2]-=region[2];
					else if(dr[2]<-regionH[2])
						dr[2]+=region[2];
					
					dr[3]=rz[j1]-rz[j2];
					if(dr[3]>regionH[3])
						dr[3]-=region[3];
					else if(dr[3]<-regionH[3])
						dr[3]+=region[3];

					g =(rvx[j1]-rvx[j2])*dr[1];
					g+=(rvy[j1]-rvy[j2])*dr[2];
					g+=(rvz[j1]-rvz[j2])*dr[3];
					g/=(2.0f*ShakeLength[nn][i]);
					ga=(float)fabs(g);
					if(ga>cDev)
						cDev=ga;
					if(ga>constraintPrec)
					{
						changed=1;
						rvx[j1]-=g*dr[1];
						rvx[j2]+=g*dr[1];
						rvy[j1]-=g*dr[2];
						rvy[j2]+=g*dr[2];
						rvz[j1]-=g*dr[3];
						rvz[j2]+=g*dr[4];
					}
				}
			}
			if(cDev>cDevV)
				cDevV=cDev;

			if(cycleR>=maxCycle) 
				printf("shake V mistake !! \n");
		}
	}
	return 0;
}


int ComputeLangevin()
{
	int n;
	for(n=1;n<=nAtom;n++)
	{
		rax[n]-=rvx[n]*Langevin;
		rax[n]+=Langevin*(float)sqrt(6*coup_temperature/deltaT)*(2.0f*(float)ran(&seed)-1.0f);
		ray[n]-=rvy[n]*Langevin;
		ray[n]+=Langevin*(float)sqrt(6*coup_temperature/deltaT)*(2.0f*(float)ran(&seed)-1.0f);
		raz[n]-=rvz[n]*Langevin;
		raz[n]+=Langevin*(float)sqrt(6*coup_temperature/deltaT)*(2.0f*(float)ran(&seed)-1.0f);
	}
	return 0;
}


int ComputeForce_hardcore() 
{
	float dr[4], tempr1,tempr2,tempr3;
	float f[4],fcVal,rr;
	float rri,rri3;
	int   j1,j2,n;	

#ifdef HBOND
	int   j0,j3;
	float rri2,rri6;
	float cos1,cos2,coscos;
	float dr0[4],dr1[4],dr2[4];
	float dr0_length,dr1_length,dr2_length;
	float cos1cos1,cos2cos2;
	float uuu,aa,a11,a22,a12;
	float f1[4],f2[4];
#endif
	
	for(j1=nstart;j1<=nend;j1++)
	{
		tempr1=rx[j1];
		tempr2=ry[j1];
		tempr3=rz[j1];

		for(n=nebrTab1[j1];n<nebrTab1[j1+1];n++)
		{
			j2=nebrTab2[n];	
			dr[1]=tempr1-rx[j2];				
			dr[2]=tempr2-ry[j2];			
			dr[3]=tempr3-rz[j2];			
			rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];
			if(rr<rrCut)
			{	
#ifdef HARDCORE
				rri=1.0f/rr;
				rri3=rri*rri*rri;
				fcVal=12*rri3*rri3*rri;
	
				fcVal*=vmForceMatrix[type[j1]][type[j2]];
				fcVal*=vmDeltaMatrix[type[j1]][type[j2]];
				
				f[1]=fcVal*dr[1];
				f[2]=fcVal*dr[2];
				f[3]=fcVal*dr[3];

				uVal+=fcVal/12*rr;

				rax[j1]+=f[1];
				rax[j2]-=f[1];
				ray[j1]+=f[2];
				ray[j2]-=f[2];
				raz[j1]+=f[3];
				raz[j2]-=f[3];
#endif
				
#ifdef HBOND
				//calculate H-bond force !!)
				if((type[j1]==1 && type[j2]==2) || (type[j1]==2 && type[j2]==1) )
				{
					j0=j1+1;  //j0: N, j1: H, j2: O, j4: 3 or exchange H and O
					j3=j2+1;


					dr0[1]=rx[j1]-rx[j0];
					if(dr0[1]>regionH[1])
						dr0[1]-=region[1];
					else if(dr0[1]<-regionH[1])
						dr0[1]+=region[1];
							
					dr0[2]=ry[j1]-ry[j0];
					if(dr0[2]>regionH[2])
						dr0[2]-=region[2];
					else if(dr0[2]<-regionH[2])
						dr0[2]+=region[2];
							
					dr0[3]=rz[j1]-rz[j0];
					if(dr0[3]>regionH[3])
						dr0[3]-=region[3];
					else if(dr0[3]<-regionH[3])
						dr0[3]+=region[3]; 							
							
					
					dr1[1]=-dr[1];
					dr1[2]=-dr[2];
					dr1[3]=-dr[3];

					
					dr2[1]=rx[j3]-rx[j2];
					if(dr2[1]>regionH[1])
						dr2[1]-=region[1];
					else if(dr2[1]<-regionH[1])
						dr2[1]+=region[1];
						
					dr2[2]=ry[j3]-ry[j2];
					if(dr2[2]>regionH[2])
						dr2[2]-=region[2];
					else if(dr2[2]<-regionH[2])
						dr2[2]+=region[2];
							
					dr2[3]=rz[j3]-rz[j2];
					if(dr2[3]>regionH[3])
						dr2[3]-=region[3];
					else if(dr2[3]<-regionH[3])
						dr2[3]+=region[3];

					cos1=0.0f;
					cos1=-(dr0[1]*dr1[1]+dr0[2]*dr1[2]+dr0[3]*dr1[3]);

					if(cos1<0)
					{
						cos2=-(dr2[1]*dr1[1]+dr2[2]*dr1[2]+dr2[3]*dr1[3]);	
						if(cos2<0)
						{	
							//calculate distance force
							dr0_length=(float)sqrt(dr0[1]*dr0[1]+dr0[2]*dr0[2]+dr0[3]*dr0[3]);
							dr1_length=(float)sqrt(rr);
							dr2_length=(float)sqrt(dr2[1]*dr2[1]+dr2[2]*dr2[2]+dr2[3]*dr2[3]);

							cos1/=dr0_length*dr1_length;
							cos2/=dr2_length*dr1_length;

							if(cos1>1.0) cos1=1.0f;
							else if(cos1<-1.0) cos1=-1.0f;
							if(cos2>1.0) cos2=1.0f;
							else if(cos2<-1.0) cos2=-1.0f;
							
							cos1cos1=cos1*cos1;
							cos2cos2=cos2*cos2;

							rri=1.0f/rr;
							rri2=rri*rri;
							rri6=rri2*rri2*rri2;
							fcVal=60*rri6*(rri-1)*H_Force;
							coscos=cos1cos1*cos2cos2;
							fcVal*=coscos;
							uuu=rri6*(5-6*rr)*H_Force;

							f[1]=fcVal*dr1[1];
							f[2]=fcVal*dr1[2];
							f[3]=fcVal*dr1[3];

							rax[j1]-=f[1];
							rax[j2]+=f[1];
							ray[j1]-=f[2];
							ray[j2]+=f[2];
							raz[j1]-=f[3];
							raz[j2]+=f[3];

							//calculate N-H-O angle potential!
							aa=uuu*cos1cos1*2*cos2;
							a11=aa*cos1/(dr0_length*dr0_length);
							a12=-aa/(dr0_length*dr1_length);
							a22=aa*cos1/(dr1_length*dr1_length);

							f1[1]=-a11*dr0[1]+a12*dr1[1];
							f2[1]=a22*dr1[1]-a12*dr0[1];
							f1[2]=-a11*dr0[2]+a12*dr1[2];
							f2[2]=a22*dr1[2]-a12*dr0[2];
							f1[3]=-a11*dr0[3]+a12*dr1[3];
							f2[3]=a22*dr1[3]-a12*dr0[3];

							rax[j0]+=f1[1];
							rax[j2]+=f2[1];
							rax[j1]-=(f1[1]+f2[1]);
							ray[j0]+=f1[2];
							ray[j2]+=f2[2];
							ray[j1]-=(f1[2]+f2[2]);
							raz[j0]+=f1[3];
							raz[j2]+=f2[3];
							raz[j1]-=(f1[3]+f2[3]);

							//calculate H-O-C' angle potential!
							aa=uuu*cos2cos2*2*cos1;
							a11=aa*cos2/(dr1_length*dr1_length);
							a12=-aa/(dr1_length*dr2_length);
							a22=aa*cos2/(dr2_length*dr2_length);

							f1[1]=-a11*dr1[1]+a12*dr2[1];
							f2[1]=a22*dr2[1]-a12*dr1[1];
							f1[2]=-a11*dr1[2]+a12*dr2[2];
							f2[2]=a22*dr2[2]-a12*dr1[2];
							f1[3]=-a11*dr1[3]+a12*dr2[3];
							f2[3]=a22*dr2[3]-a12*dr1[3];

							rax[j1]+=f1[1];
							rax[j3]+=f2[1];
							rax[j2]-=(f1[1]+f2[1]);
							ray[j1]+=f1[2];
							ray[j3]+=f2[2];
							ray[j2]-=(f1[2]+f2[2]);
							raz[j1]+=f1[3];
							raz[j3]+=f2[3];
							raz[j2]-=(f1[3]+f2[3]);	 

							//calculate H-bond energy and numebr
							hVal+=uuu*coscos;
							hnumberVal+=1;	
						}
					}
				}// H-BOND
#endif
			}//if rr<rrCut
		 } //for(n=nebrTab1[j1];n<nebrTab1[j1+1];n++)
		
		for(n=nebrTab3[j1];n<nebrTab3[j1+1];n++)
		{
			j2=nebrTab4[n];
			
			dr[1]=tempr1-rx[j2];
			if(fabs(dr[1])>regionH[1])
			{
				if(dr[1]>0)
					dr[1]-=region[1];
				else 
					dr[1]+=region[1];
			}
				
			dr[2]=tempr2-ry[j2];
			if(fabs(dr[2])>regionH[2])
			{
				if(dr[2]>0)
					dr[2]-=region[2];
				else 
					dr[2]+=region[2];
			}	
			
			dr[3]=tempr3-rz[j2];			
			if(fabs(dr[3])>regionH[3])
			{
				if(dr[3]>0)
					dr[3]-=region[3];
				else 
					dr[3]+=region[3];
			}
			
			rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];	

			if(rr<rrCut)
			{
#ifdef HARDCORE
				rri=1.0f/rr;
				rri3=rri*rri*rri;
				fcVal=12*rri3*rri3*rri;
	
				fcVal*=vmForceMatrix[type[j1]][type[j2]];
				fcVal*=vmDeltaMatrix[type[j1]][type[j2]];
				
				f[1]=fcVal*dr[1];
				f[2]=fcVal*dr[2];
				f[3]=fcVal*dr[3];

				uVal+=fcVal/12*rr;

				rax[j1]+=f[1];
				rax[j2]-=f[1];
				ray[j1]+=f[2];
				ray[j2]-=f[2];
				raz[j1]+=f[3];
				raz[j2]-=f[3];
#endif
#ifdef HBOND				  
				//calculate H-bond force !!)
				if((type[j1]==1 && type[j2]==2) || (type[j1]==2 && type[j2]==1) )
				{
					j0=j1+1;  //j0: N, j1: H, j2: O, j4: 3 or exchange H and O
					j3=j2+1;


					dr0[1]=rx[j1]-rx[j0];
					if(dr0[1]>regionH[1])
						dr0[1]-=region[1];
					else if(dr0[1]<-regionH[1])
						dr0[1]+=region[1];
							
					dr0[2]=ry[j1]-ry[j0];
					if(dr0[2]>regionH[2])
						dr0[2]-=region[2];
					else if(dr0[2]<-regionH[2])
						dr0[2]+=region[2];
							
					dr0[3]=rz[j1]-rz[j0];
					if(dr0[3]>regionH[3])
						dr0[3]-=region[3];
					else if(dr0[3]<-regionH[3])
						dr0[3]+=region[3]; 							
							
					
					dr1[1]=-dr[1];
					dr1[2]=-dr[2];
					dr1[3]=-dr[3];

					
					dr2[1]=rx[j3]-rx[j2];
					if(dr2[1]>regionH[1])
						dr2[1]-=region[1];
					else if(dr2[1]<-regionH[1])
						dr2[1]+=region[1];
						
					dr2[2]=ry[j3]-ry[j2];
					if(dr2[2]>regionH[2])
						dr2[2]-=region[2];
					else if(dr2[2]<-regionH[2])
						dr2[2]+=region[2];
							
					dr2[3]=rz[j3]-rz[j2];
					if(dr2[3]>regionH[3])
						dr2[3]-=region[3];
					else if(dr2[3]<-regionH[3])
						dr2[3]+=region[3];

					cos1=0.0f;
					cos1=-(dr0[1]*dr1[1]+dr0[2]*dr1[2]+dr0[3]*dr1[3]);

					if(cos1<0)
					{
						cos2=-(dr2[1]*dr1[1]+dr2[2]*dr1[2]+dr2[3]*dr1[3]);	
						if(cos2<0)
						{	 
							//calculate distance force
							dr0_length=(float)sqrt(dr0[1]*dr0[1]+dr0[2]*dr0[2]+dr0[3]*dr0[3]);
							dr1_length=(float)sqrt(rr);
							dr2_length=(float)sqrt(dr2[1]*dr2[1]+dr2[2]*dr2[2]+dr2[3]*dr2[3]);

							cos1/=dr0_length*dr1_length;
							cos2/=dr2_length*dr1_length;

							if(cos1>1.0) cos1=1.0f;
							else if(cos1<-1.0) cos1=-1.0f;
							if(cos2>1.0) cos2=1.0f;
							else if(cos2<-1.0) cos2=-1.0f;
							
							cos1cos1=cos1*cos1;
							cos2cos2=cos2*cos2;

							
							rri=1.0f/rr;
							rri2=rri*rri;
							rri6=rri2*rri2*rri2;
							fcVal=60*rri6*(rri-1)*H_Force;							
							coscos=cos1cos1*cos2cos2;
							fcVal*=coscos;
							uuu=rri6*(5-6*rr)*H_Force;

							f[1]=fcVal*dr1[1];
							f[2]=fcVal*dr1[2];
							f[3]=fcVal*dr1[3];

							rax[j1]-=f[1];
							rax[j2]+=f[1];
							ray[j1]-=f[2];
							ray[j2]+=f[2];
							raz[j1]-=f[3];
							raz[j2]+=f[3];

							//calculate N-H-O angle potential!
							aa=uuu*cos1cos1*2*cos2;
							a11=aa*cos1/(dr0_length*dr0_length);
							a12=-aa/(dr0_length*dr1_length);
							a22=aa*cos1/(dr1_length*dr1_length);

							f1[1]=-a11*dr0[1]+a12*dr1[1];
							f2[1]=a22*dr1[1]-a12*dr0[1];
							f1[2]=-a11*dr0[2]+a12*dr1[2];
							f2[2]=a22*dr1[2]-a12*dr0[2];
							f1[3]=-a11*dr0[3]+a12*dr1[3];
							f2[3]=a22*dr1[3]-a12*dr0[3];

							rax[j0]+=f1[1];
							rax[j2]+=f2[1];
							rax[j1]-=(f1[1]+f2[1]);
							ray[j0]+=f1[2];
							ray[j2]+=f2[2];
							ray[j1]-=(f1[2]+f2[2]);
							raz[j0]+=f1[3];
							raz[j2]+=f2[3];
							raz[j1]-=(f1[3]+f2[3]);

							//calculate H-O-C' angle potential!
							aa=uuu*cos2cos2*2*cos1;
							a11=aa*cos2/(dr1_length*dr1_length);
							a12=-aa/(dr1_length*dr2_length);
							a22=aa*cos2/(dr2_length*dr2_length);

							f1[1]=-a11*dr1[1]+a12*dr2[1];
							f2[1]=a22*dr2[1]-a12*dr1[1];
							f1[2]=-a11*dr1[2]+a12*dr2[2];
							f2[2]=a22*dr2[2]-a12*dr1[2];
							f1[3]=-a11*dr1[3]+a12*dr2[3];
							f2[3]=a22*dr2[3]-a12*dr1[3];

							rax[j1]+=f1[1];
							rax[j3]+=f2[1];
							rax[j2]-=(f1[1]+f2[1]);
							ray[j1]+=f1[2];
							ray[j3]+=f2[2];
							ray[j2]-=(f1[2]+f2[2]);
							raz[j1]+=f1[3];
							raz[j3]+=f2[3];
							raz[j2]-=(f1[3]+f2[3]);	 

							//calculate H-bond energy and numebr
							hVal+=uuu*coscos;
							hnumberVal+=1;

						}
					}
				}//h-bond
#endif

			}//if rr<rrCut
		 } //for(n=nebrTab1[j1];n<nebrTab1[j1+1];n++)
	}//for (j1=

	return 0;
}

int   ComputeForce_hydrophobic()
{
	float dr[4], tempr1,tempr2,tempr3;
	float f[4],fcVal,rr;
	float rri,rri3;
	int   j1,j2,n; 


	for(j1=nstart;j1<=nend;j1++)
	{
		if(type[j1]!=5)
			continue;
		tempr1=rx[j1];
		tempr2=ry[j1];
		tempr3=rz[j1];

		for(n=nebrTab1_second[j1];n<nebrTab1_second[j1+1];n++)
		{
			j2=nebrTab2_second[n];
			if(type2[j1]==1	&& type2[j2] ==1 )
			{
			
				dr[1]=tempr1-rx[j2];				
				dr[2]=tempr2-ry[j2];			
				dr[3]=tempr3-rz[j2];			
				rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];
				if(rr<rrCut_second)
				{	 
					rri=1.0f/rr;
					rri3=rri*rri*rri;

					fcVal=12*rri3*rri*244.140625f*(rri3*244.140625f-1.0f);
					fcVal*=Hydro_Force;

					hydroVal+=Hydro_Force*(rri3*59604.644775390625f-2*244.140625f)*rri3;

					f[1]=fcVal*dr[1];
					f[2]=fcVal*dr[2];
					f[3]=fcVal*dr[3];

					rax[j1]+=f[1];
					rax[j2]-=f[1];
					ray[j1]+=f[2];
					ray[j2]-=f[2];
					raz[j1]+=f[3];
					raz[j2]-=f[3];					

				}//if rr<rrCut
			}
		 } //for(n=nebrTab1[j1];n<nebrTab1[j1+1];n++)
		
		for(n=nebrTab3_second[j1];n<nebrTab3_second[j1+1];n++)
		{
			j2=nebrTab4_second[n];
 			if(type2[j1]==1	&& type2[j2] ==1 )
			{

				dr[1]=tempr1-rx[j2];
				if(fabs(dr[1])>regionH[1])
				{
					if(dr[1]>0)
						dr[1]-=region[1];
					else 
						dr[1]+=region[1];
				}
					
				dr[2]=tempr2-ry[j2];
				if(fabs(dr[2])>regionH[2])
				{
					if(dr[2]>0)
						dr[2]-=region[2];
					else 
						dr[2]+=region[2];
				}	
				
				dr[3]=tempr3-rz[j2];			
				if(fabs(dr[3])>regionH[3])
				{
					if(dr[3]>0)
						dr[3]-=region[3];
					else 
						dr[3]+=region[3];
				}
				
				rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];	

				if(rr<rrCut_second)
				{		 
					rri=1.0f/rr;
					rri3=rri*rri*rri;

					fcVal=12*rri3*rri*(rri3*59604.644775390625f-244.140625f);
					fcVal*=Hydro_Force;

					hydroVal+=Hydro_Force*(rri3*59604.644775390625f-2*244.140625f)*rri3;

					f[1]=fcVal*dr[1];
					f[2]=fcVal*dr[2];
					f[3]=fcVal*dr[3];

					rax[j1]+=f[1];
					rax[j2]-=f[1];
					ray[j1]+=f[2];
					ray[j2]-=f[2];
					raz[j1]+=f[3];
					raz[j2]-=f[3];
				}//if rr<rrCut
			}
		 } //for(n=nebrTab1[j1];n<nebrTab1[j1+1];n++)
	}//for (j1=	
	return 0;
}



int  ComputeBondForce()  //compute bond length and bond angle force!!
{
	float dr[3][4],f[4];
	int   i,j1,j2,n,nn,k;

#ifdef BOND
	float r2,r1;
#endif

#ifdef BONDANGLE
	float ab,aa;
	float f1[4],f2[4];
	float tempdr[4];
	float tempf1[4],tempf2[4],tempf3[4];
	int j3;
	float radii1,radii2,a11,a12,a22;
	float rr[3];
#endif

#ifdef TORSION
	int j4;
	float sb1,sb2,sb3,rb1,rb3;
	float a,c,c0,b1mag,b2mag,b3mag,b1mag2,b2mag2,b3mag2;
	float ctmp,r12c1,r12c2,c1mag,c2mag;
	float sc1,sc2,s1,s2,s12;
	float sx1,sx2,sx12,sy1,sy2,sy12,sz1,sz2,sz12;
	float a13,a23,a33;
	float cc;
#endif

//1-4 interaction
#ifdef INTERACTION14
	float dr14[4];  
	float rr14;
	float rri,rri3;
	float fcVal;
#endif

#ifdef OUTOFPLANE
	float f3[4];
	float p;
#endif


	for(nn=1;nn<=TypePoly;nn++)
	{

#ifdef MYMPI
		for(n=nchainstart[nn];n<=nchainend[nn];n++)
		{

#else
		for(n=1;n<=nChain[nn];n++)
		{
#endif

#ifdef BOND
			for(i=1;i<=BondNumber[nn];i++)
			{
				j1=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondIndex1[nn][i];
				j2=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondIndex2[nn][i];
				
				dr[1][1]=rx[j1]-rx[j2];
				if(dr[1][1]>regionH[1])
					dr[1][1]-=region[1];
				else if(dr[1][1]<-regionH[1])
					dr[1][1]+=region[1];
				
				dr[1][2]=ry[j1]-ry[j2];
				if(dr[1][2]>regionH[2])
					dr[1][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[1][2]+=region[2];
				
				dr[1][3]=rz[j1]-rz[j2];
				if(dr[1][3]>regionH[3])
					dr[1][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[1][3]+=region[3];

				r1=dr[1][1]*dr[1][1]+dr[1][2]*dr[1][2]+dr[1][3]*dr[1][3];
				r2=(float)(sqrt(r1));
#ifdef	MYDEBUG
				if(r2>regionH[2])
				{
					printf("Warning, bond length too long, %d and %d is %f at step %d \r\n",j1,j2,r2,stepCount); 					
				}

				else if(r2<0.04)
				{
					printf("Warning, bond length too short, %d and %d is %f at step %d \r\n",
						j1,j2,r2,stepCount);
				}
#endif
				for(k=1;k<=3;k++)
				{
					f[k]=(float)(BondStrenth[nn][i]*(BondLength[nn][i]/r2-1)*dr[1][k]);
				}
				rax[j1]+=f[1];
				rax[j2]-=f[1];
				ray[j1]+=f[2];
				ray[j2]-=f[2];
				raz[j1]+=f[3];
				raz[j2]-=f[3];

#ifdef PRESSURE
				int kk;
				for(kk=1;kk<=3;kk++)
				for(k=1;k<=3;k++)
				{
					virSum[k][kk]+=f[k]*dr[1][kk];
				}
#endif

#ifdef GETTENSION
				Gettension(2,f,dr[1],j1,j2);
#endif
				virSum[0][0]+=BondStrenth[nn][i]*(BondLength[nn][i]-r2)*r2;
				bondVal+=(float)(BondStrenth[nn][i]/2.0*(r2-BondLength[nn][i])*(r2-BondLength[nn][i]));
			}

#endif


#ifdef BONDANGLE
			for(i=1;i<=BondAngleNumber[nn];i++)
			{
				j1=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondAngleIndex1[nn][i];
				j2=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondAngleIndex2[nn][i];
				j3=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondAngleIndex3[nn][i];

				dr[1][1]=rx[j1]-rx[j2];
				if(dr[1][1]>regionH[1])
					dr[1][1]-=region[1];
				else if(dr[1][1]<-regionH[1])
					dr[1][1]+=region[1];
				
				dr[1][2]=ry[j1]-ry[j2];
				if(dr[1][2]>regionH[2])
					dr[1][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[1][2]+=region[2];
				
				dr[1][3]=rz[j1]-rz[j2];
				if(dr[1][3]>regionH[3])
					dr[1][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[1][3]+=region[3];
				
				
				dr[2][1]=rx[j3]-rx[j2];
				if(dr[2][1]>regionH[1])
					dr[2][1]-=region[1];
				else if(dr[2][1]<-regionH[1])
					dr[2][1]+=region[1];
				
				dr[2][2]=ry[j3]-ry[j2];
				if(dr[2][2]>regionH[2])
					dr[2][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[2][2]+=region[2];
				
				dr[2][3]=rz[j3]-rz[j2];
				if(dr[2][3]>regionH[3])
					dr[2][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[2][3]+=region[3];								


				
				rr[1]=dr[1][1]*dr[1][1]+dr[1][2]*dr[1][2]+dr[1][3]*dr[1][3];
				rr[2]=dr[2][1]*dr[2][1]+dr[2][2]*dr[2][2]+dr[2][3]*dr[2][3];
				radii1=(float)sqrt(rr[1]);
				radii2=(float)sqrt(rr[2]);

				ab= dr[1][1]*dr[2][1]+dr[1][2]*dr[2][2]+dr[1][3]*dr[2][3];
				ab/=radii1*radii2;
				
				if(ab>1.0) ab=1.0f;
				else if(ab<-1.0) ab=-1.0f;			

				aa=ab-(float)(cos(BondAngle[nn][i]));
				a11=aa*BondAngleStrenth[nn][i]*ab/rr[1];
				a12=-aa*BondAngleStrenth[nn][i]/(radii1*radii2);
				a22=aa*BondAngleStrenth[nn][i]*ab/rr[2];

				f1[1]=a11*dr[1][1]+a12*dr[2][1];
				f2[1]=a22*dr[2][1]+a12*dr[1][1];
				f1[2]=a11*dr[1][2]+a12*dr[2][2];
				f2[2]=a22*dr[2][2]+a12*dr[1][2];
				f1[3]=a11*dr[1][3]+a12*dr[2][3];
				f2[3]=a22*dr[2][3]+a12*dr[1][3];

				rax[j1]+=f1[1];
				rax[j3]+=f2[1];
				rax[j2]-=(f1[1]+f2[1]);
				ray[j1]+=f1[2];
				ray[j3]+=f2[2];
				ray[j2]-=(f1[2]+f2[2]);
				raz[j1]+=f1[3];
				raz[j3]+=f2[3];
				raz[j2]-=(f1[3]+f2[3]);

				for(k=1;k<=3;k++)
				{
					
					tempdr[k]=dr[1][k]-dr[2][k];
					tempf1[k]=(f1[k]*2+f2[k])/3;
					tempf2[k]=(f2[k]*2+f1[k])/3;
					tempf3[k]=(f1[k]-f2[k])/3;

					virSum[0][0]+=tempf1[k]*dr[1][k];
					virSum[0][0]+=tempf2[k]*dr[2][k];
					virSum[0][0]+=tempf3[k]*tempdr[k];
				}
#ifdef PRESSURE
				int kk;
				for(kk=1;kk<=3;kk++)
				for(k=1;k<=3;k++)
				{
					virSum[k][kk]+=tempf1[k]*dr[1][kk];
					virSum[k][kk]+=tempf2[k]*dr[2][kk];
					virSum[k][kk]+=tempf3[k]*tempdr[kk];
				}
#endif

#ifdef GETTENSION
				Gettension(3,tempf1, dr[1],j1,j2);
				Gettension(3,tempf2, dr[2],j2,j3);
				Gettension(3,tempf3,tempdr,j1,j3);
#endif
			}
#endif //bond angle

#ifdef TORSION
			for(i=1;i<=BondTorsionNumber[nn];i++)
			{
				j1=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondTorsionIndex1[nn][i];
				j2=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondTorsionIndex2[nn][i];
				j3=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondTorsionIndex3[nn][i];
				j4=ChainBegin[nn]+(n-1)*ChainLength[nn]+BondTorsionIndex4[nn][i];
			
				dr[1][1]=rx[j1]-rx[j2];
				if(dr[1][1]>regionH[1])
					dr[1][1]-=region[1];
				else if(dr[1][1]<-regionH[1])
					dr[1][1]+=region[1];
				
				dr[1][2]=ry[j1]-ry[j2];
				if(dr[1][2]>regionH[2])
					dr[1][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[1][2]+=region[2];
				
				dr[1][3]=rz[j1]-rz[j2];
				if(dr[1][3]>regionH[3])
					dr[1][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[1][3]+=region[3];
				
				
				dr[2][1]=rx[j3]-rx[j2];
				if(dr[2][1]>regionH[1])
					dr[2][1]-=region[1];
				else if(dr[2][1]<-regionH[1])
					dr[2][1]+=region[1];
				
				dr[2][2]=ry[j3]-ry[j2];
				if(dr[2][2]>regionH[2])
					dr[2][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[2][2]+=region[2];
				
				dr[2][3]=rz[j3]-rz[j2];
				if(dr[2][3]>regionH[3])
					dr[2][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[2][3]+=region[3];
				
				
				dr[3][1]=rx[j4]-rx[j3];
				if(dr[3][1]>regionH[1])
					dr[3][1]-=region[1];
				else if(dr[3][1]<-regionH[1])
					dr[3][1]+=region[1];
				
				dr[3][2]=ry[j4]-ry[j3];
				if(dr[3][2]>regionH[2])
					dr[3][2]-=region[2];
				else if(dr[3][2]<-regionH[2])
					dr[3][2]+=region[2];
				
				dr[3][3]=rz[j4]-rz[j3];
				if(dr[3][3]>regionH[3])
					dr[3][3]-=region[3];
				else if(dr[3][3]<-regionH[3])
					dr[3][3]+=region[3];

				sb1 = 1.0f / (dr[1][1]*dr[1][1] + dr[1][2]*dr[1][2] + dr[1][3]*dr[1][3]);
				sb2 = 1.0f / (dr[2][1]*dr[2][1] + dr[2][2]*dr[2][2] + dr[2][3]*dr[2][3]);
				sb3 = 1.0f / (dr[3][1]*dr[3][1] + dr[3][2]*dr[3][2] + dr[3][3]*dr[3][3]);
        
				rb1 = (float)sqrt(sb1);
				rb3 = (float)sqrt(sb3);
        
				c0 = (dr[1][1]*dr[3][1] + dr[1][2]*dr[3][2] + dr[1][3]*dr[3][3]) * rb1*rb3;

				// 1st and 2nd angle
        
				b1mag2 = dr[1][1]*dr[1][1] + dr[1][2]*dr[1][2] + dr[1][3]*dr[1][3];
				b1mag = (float)sqrt(b1mag2);
				b2mag2 = dr[2][1]*dr[2][1] + dr[2][2]*dr[2][2] + dr[2][3]*dr[2][3];
				b2mag = (float)sqrt(b2mag2);
				b3mag2 = dr[3][1]*dr[3][1] + dr[3][2]*dr[3][2] + dr[3][3]*dr[3][3];
				b3mag = (float)sqrt(b3mag2);

				ctmp = dr[1][1]*dr[2][1] + dr[1][2]*dr[2][2] + dr[1][3]*dr[2][3];
				r12c1 = 1.0f / (b1mag*b2mag);
				c1mag = ctmp * r12c1;

				ctmp = dr[3][1]*dr[2][1] + dr[3][2]*dr[2][2] + dr[3][3]*dr[2][3];
				r12c2 = 1.0f / (b2mag*b3mag);
				c2mag = ctmp * r12c2;

				// cos and sin of 2 angles and final c

				sc1 = (float)sqrt(1.0 - c1mag*c1mag);
				if (sc1 < SMALL) sc1 = SMALL;
				sc1 = 1.0f/sc1;

				sc2 = (float)sqrt(1.0 - c2mag*c2mag);
				if (sc2 < SMALL) sc2 = SMALL;
				sc2 = 1.0f/sc2;

				s1 = sc1 * sc1;
				s2 = sc2 * sc2;
				s12 = sc1 * sc2;
				c = (c0 - c1mag*c2mag) * s12;
    
				if (c > 1.0f) 
					c = 1.0f;
				if (c < -1.0f)
					c = -1.0f;
			
				 //Calculate phi and psi distribution here!!
				/*
				a=(dr[2][2]*dr[1][3]-dr[2][3]*dr[1][2])*dr[3][1]
					+ (dr[2][3]*dr[1][1]-dr[2][1]*dr[1][3])*dr[3][2]
					+(dr[2][1]*dr[1][2]-dr[2][2]*dr[1][1])*dr[3][3];
				   				
				fa=(float)acos(c);  // fa is between 0 and Pi
				if(a<0) fa=-fa;

				if(i!=1 && i%3==1)	 //psi
				{
					fout_dis<<" "<<fa*180/3.14159f<<"\n";
				}				
				else if(i!=45 && i%3==0)	//phi
				{
					fout_dis<<fa*180/3.14159f;;
				}*/ 

				a = BondTorsionStrenth[nn][i];
				cc=c*c; 

				if(BondTorsionType[nn][i]==1)  //u=cos(3pha);
				{
					a*=1.5f*(4*cc-1);
					a*=(float)sqrt(1-cc);
				}

				else if(BondTorsionType[nn][i]==2)  //u=cos( 3(fa-120) ); 
				{
					a*=1.5f*(4*cc-1);
					a*=(float)sqrt(1-cc);
				} 
				else //u=cos(fa)
				{
					a*=(float)sqrt(1-cc);
					a/=2;
				}

			    c = c * a;
				s12 = s12 * a;

				a11 = (-c*sb1*s1);
				a22 = sb2*(2.0f*c0*s12 - c*(s1+s2));
				a33 = (-c*sb3*s2);
				a12 = r12c1*(c1mag*c*s1 + c2mag*s12);
				a13 = rb1*rb3*s12;
				a23 = r12c2*(-c2mag*c*s2 - c1mag*s12);

				sx1  = a11*dr[1][1] + a12*dr[2][1] + a13*dr[3][1];
				sx2  = a12*dr[1][1] + a22*dr[2][1] + a23*dr[3][1];
				sx12 = a13*dr[1][1] + a23*dr[2][1] + a33*dr[3][1];
				sy1  = a11*dr[1][2] + a12*dr[2][2] + a13*dr[3][2];
				sy2  = a12*dr[1][2] + a22*dr[2][2] + a23*dr[3][2];
				sy12 = a13*dr[1][2] + a23*dr[2][2] + a33*dr[3][2];
				sz1  = a11*dr[1][3] + a12*dr[2][3] + a13*dr[3][3];
				sz2  = a12*dr[1][3] + a22*dr[2][3] + a23*dr[3][3];
				sz12 = a13*dr[1][3] + a23*dr[2][3] + a33*dr[3][3];

				rax[j1] -= sx1;
				ray[j1] -= sy1;
				raz[j1] -= sz1;
				
				rax[j2] += (sx2 + sx1);
				ray[j2] += (sy2 + sy1);
				raz[j2] += (sz2 + sz1);
				
				rax[j3] += (sx12 - sx2);
				ray[j3] += (sy12 - sy2);
				raz[j3] += (sz12 - sz2);
				
				rax[j4] -= sx12;
				ray[j4] -= sy12;
				raz[j4] -= sz12;
			}

#endif  //torsion
#ifdef  OUTOFPLANE
			for(i=1;i<=OutofPlaneNumber[nn];i++)
			{
			   	j1=ChainBegin[nn]+(n-1)*ChainLength[nn]+OutofPlaneIndex1[nn][i];
				j2=ChainBegin[nn]+(n-1)*ChainLength[nn]+OutofPlaneIndex2[nn][i];
				j3=ChainBegin[nn]+(n-1)*ChainLength[nn]+OutofPlaneIndex3[nn][i];
				j4=ChainBegin[nn]+(n-1)*ChainLength[nn]+OutofPlaneIndex4[nn][i];

				dr[1][1]=rx[j1]-rx[j4];
				if(dr[1][1]>regionH[1])
					dr[1][1]-=region[1];
				else if(dr[1][1]<-regionH[1])
					dr[1][1]+=region[1];
				
				dr[1][2]=ry[j1]-ry[j4];
				if(dr[1][2]>regionH[2])
					dr[1][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[1][2]+=region[2];
				
				dr[1][3]=rz[j1]-rz[j4];
				if(dr[1][3]>regionH[3])
					dr[1][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[1][3]+=region[3];
				
				
				dr[2][1]=rx[j2]-rx[j4];
				if(dr[2][1]>regionH[1])
					dr[2][1]-=region[1];
				else if(dr[2][1]<-regionH[1])
					dr[2][1]+=region[1];
				
				dr[2][2]=ry[j2]-ry[j4];
				if(dr[2][2]>regionH[2])
					dr[2][2]-=region[2];
				else if(dr[1][2]<-regionH[2])
					dr[2][2]+=region[2];
				
				dr[2][3]=rz[j2]-rz[j4];
				if(dr[2][3]>regionH[3])
					dr[2][3]-=region[3];
				else if(dr[1][3]<-regionH[3])
					dr[2][3]+=region[3];
				
				
				dr[3][1]=rx[j3]-rx[j4];
				if(dr[3][1]>regionH[1])
					dr[3][1]-=region[1];
				else if(dr[3][1]<-regionH[1])
					dr[3][1]+=region[1];
				
				dr[3][2]=ry[j3]-ry[j4];
				if(dr[3][2]>regionH[2])
					dr[3][2]-=region[2];
				else if(dr[3][2]<-regionH[2])
					dr[3][2]+=region[2];
				
				dr[3][3]=rz[j3]-rz[j4];
				if(dr[3][3]>regionH[3])
					dr[3][3]-=region[3];
				else if(dr[3][3]<-regionH[3])
					dr[3][3]+=region[3];
				
				p=(dr[1][2]*dr[2][3]-dr[1][3]*dr[2][2])*dr[3][1]
					+(dr[1][3]*dr[2][1]-dr[1][1]*dr[2][3])*dr[3][2]
					+(dr[1][1]*dr[2][2]-dr[1][2]*dr[2][1])*dr[3][3];

				if(p<0)
				{
					outofplaneNumber++;

					f3[1]=(dr[1][2]*dr[2][3]-dr[1][3]*dr[2][2])*OutofPlaneStrenth[n][i];
					f3[2]=(dr[1][3]*dr[2][1]-dr[1][1]*dr[2][3])*OutofPlaneStrenth[n][i];
					f3[3]=(dr[1][1]*dr[2][2]-dr[1][2]*dr[2][1])*OutofPlaneStrenth[n][i];
					
					f2[1]=(dr[3][2]*dr[1][3]-dr[3][3]*dr[1][2])*OutofPlaneStrenth[n][i];
					f2[2]=(dr[3][3]*dr[1][1]-dr[3][1]*dr[1][3])*OutofPlaneStrenth[n][i];
					f2[3]=(dr[3][1]*dr[1][2]-dr[3][2]*dr[1][1])*OutofPlaneStrenth[n][i];
					
					f1[1]=(dr[2][2]*dr[3][3]-dr[2][3]*dr[3][2])*OutofPlaneStrenth[n][i];
					f1[2]=(dr[2][3]*dr[3][1]-dr[2][1]*dr[3][3])*OutofPlaneStrenth[n][i];
					f1[3]=(dr[2][1]*dr[3][2]-dr[2][2]*dr[3][1])*OutofPlaneStrenth[n][i];

					rax[j1]+=f1[1];
					ray[j1]+=f1[2];
					raz[j1]+=f1[3];
					
					rax[j2]+=f2[1];
					ray[j2]+=f2[2];
					raz[j2]+=f2[3];
					
					rax[j3]+=f3[1];
					ray[j3]+=f3[2];
					raz[j3]+=f3[3];
					
					rax[j4]+=-f1[1]-f2[1]-f3[1];
					ray[j4]+=-f1[2]-f2[2]-f3[2];
					raz[j4]+=-f1[3]-f2[3]-f3[3];
				}
			}
#endif

#ifdef  INTERACTION14
			for(i=1;i<=Interaction14Number[nn];i++)
			{
				j1=ChainBegin[nn]+(n-1)*ChainLength[nn]+Interaction14Index1[nn][i];
				j2=ChainBegin[nn]+(n-1)*ChainLength[nn]+Interaction14Index2[nn][i];

				dr14[1]=rx[j1]-rx[j2];
				if(dr14[1]>regionH[1])
					dr14[1]-=region[1];
				else if(dr14[1]<-regionH[1])
					dr14[1]+=region[1];
				
				dr14[2]=ry[j1]-ry[j2];
				if(dr14[2]>regionH[2])
					dr14[2]-=region[2];
				else if(dr14[2]<-regionH[2])
					dr14[2]+=region[2];
				
				dr14[3]=rz[j1]-rz[j2];
				if(dr14[3]>regionH[3])
					dr14[3]-=region[3];
				else if(dr14[3]<-regionH[3])
					dr14[3]+=region[3];

				rr14=dr14[1]*dr14[1]+dr14[2]*dr14[2]+dr14[3]*dr14[3];	

				if(rr14<rrCut)
				{

					rri=1.0f/rr14;
					rri3=rri*rri*rri;
					fcVal=12*rri3*rri3*rri;
		
					fcVal*=vmForceMatrix14[type[j1]][type[j2]];
					fcVal*=vmDeltaMatrix14[type[j1]][type[j2]];

					f[1]=fcVal*dr14[1];
					f[2]=fcVal*dr14[2];
					f[3]=fcVal*dr14[3];

					uVal+=fcVal/12*rr14;

					rax[j1]+=f[1];
					rax[j2]-=f[1];
					ray[j1]+=f[2];
					ray[j2]-=f[2];
					raz[j1]+=f[3];
					raz[j2]-=f[3];
				}
			}// 1-4 interaction
#endif
		}// for(n=1;n<=nChain[nn];n++)
	}//for(nn=1;nn<=TypePoly;nn++)
	return 0;
}




#ifdef MYMPI
void FoldForce()
{
	MPI_Allreduce(rax+1,raxx+1,nAtom,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(ray+1,rayy+1,nAtom,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(raz+1,razz+1,nAtom,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
}

void Expandr()
{
	MPI_Status status;
	if(node==0)
	{
		MPI_Send(rx+nstart,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD);
		MPI_Recv(rx+nend+1,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD,&status );
	}
	else if(node==1)
	{	
		MPI_Recv(rx+1,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status );
		MPI_Send(rx+nstart,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
	
	if(node==0)
	{
		MPI_Send(ry+nstart,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD);
		MPI_Recv(ry+nend+1,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD,&status );
	}
	else if(node==1)
	{
		MPI_Recv(ry+1,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status );
		MPI_Send(ry+nstart,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
	
	if(node==0)
	{
		MPI_Send(rz+nstart,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD);
		MPI_Recv(rz+nend+1,nlocal,MPI_FLOAT,1,1,MPI_COMM_WORLD,&status );
	}
	else if(node==1)
	{
		MPI_Recv(rz+1,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD,&status );
		MPI_Send(rz+nstart,nlocal,MPI_FLOAT,0,1,MPI_COMM_WORLD);
	}
}

#endif


int  BuildNebrList() //nebr list method
{
#ifdef SSL
	__declspec(align(16)) float dr[4], invWid[4],shift[4];
#else
	float dr[4], invWid[4],shift[4];
#endif
	float rr,rrNebr;
	int  bTran,bTran0;
	float xtmp,ytmp,ztmp;
	int    c,j1,j2,j11,j22,k,m1X,m1Y,m1Z,m2,m2X,m2Y,m2Z,n,offset;
	int   inPolyj1;
	// below three array define the half-neighbour cell(total 14), include itself.
	int    iofX[]={0,0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1};
	int    iofY[]={0,0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
	int    iofZ[]={0,0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1};
	int    c1,c2,c3;
	int    i;

	nebrNow=0;  //clear flag for buildneighborlist
	dispHi=0.0f;   //clear the sum of max V, which is used to decide whether to buildlist or not!

#ifdef MYDEBUG
	printf("Build neighbour list when step %d, averaged number is %d \r\n", 
			stepCount,nebrTabLen/nAtom);
#endif


	rrNebr=(rCut+rNebrShell)*(rCut+rNebrShell);
	for(k=1;k<=3;k++)
		invWid[k]=cells[k]/region[k];
//below code: cellList play dual role: The first parts of the array consists of pointers linking
//different atoms belong to the same cell, while the second parts one per cell, point to the
//first atom in each cell, 0 mean the final atom in one list of cell and an empty cell respectively
	for(n=nAtom+1;n<=nAtom+cells[1]*cells[2]*cells[3];n++)
		cellList[n]=0;

	
	for(n=1;n<=nAtom;n++)
	{
		c3=(int)((rz[n]+regionH[3])*invWid[3]);
		c2=(int)((ry[n]+regionH[2])*invWid[2]);
		c1=(int)((rx[n]+regionH[1])*invWid[1]);
		if(c3==cells[3]) c3=0;
		if(c2==cells[2]) c2=0;
		if(c1==cells[1]) c1=0;
		c=((int)c3*cells[2]+(int)c2)*cells[1]+(int)c1+nAtom+1;
		if(c>nAtom+cells[1]*cells[2]*cells[3] || c<nAtom+1)
		{
			printf("In step %d for atom number %d \r\n",stepCount,n);
			printf(" 3  2  1 is %d,%d,%d\r\n", cells[3],cells[2],cells[1]);
			printf("c3,c2,c1 is %d,%d,%d\r\n",c3,c2,c1);
			printf("n is %d, rx[n] is %f, regionH[1] is %f, invWid[1] is %f\r\n",n,rx[n],regionH[1],invWid[1]);
			printf("n is %d, ry[n] is %f, regionH[2] is %f, invWid[2] is %f\r\n",n,ry[n],regionH[2],invWid[2]);
			printf("n is %d, rz[n] is %f, regionH[3] is %f, invWid[3] is %f\r\n",n,rz[n],regionH[3],invWid[3]);
			return 2;  //over flow !!!
		}
		cellList[n]=cellList[c];
		cellList[c]=n;
	}
	nebrTabLen=0;
	nebrTabLen2=0;
	
	for(j1=nstart;j1<=nend;j1++)
	{
		inPolyj1=inPoly[j1];
		if(inPolyj1!=0)
		{
			i=1;
			j11=j1;
			while(j11>nChain[i]*ChainLength[i])
			{
				i++;
				j11=j1-ChainBegin[i];								
			}
			j11=j1-ChainBegin[i];
		}

		nebrTab1[j1]=nebrTabLen+1;
		nebrTab3[j1]=nebrTabLen2+1;

		xtmp=rx[j1];
		ytmp=ry[j1];
		ztmp=rz[j1];

		bTran0=0;
		m1Z=(int)((rz[j1]+regionH[3])*invWid[3]);
		m1Y=(int)((ry[j1]+regionH[2])*invWid[2]);
		m1X=(int)((rx[j1]+regionH[1])*invWid[1]);
		if(m1Z==cells[3]) 
		{
			m1Z=0;
			bTran0=1;
		}
		else if(m1Z==cells[3]-1 || m1Z==0)
			bTran0=1;
		if(m1Y==cells[2])
		{
			m1Y=0;
			bTran0=1;
		}
		else if(m1Y==cells[2]-1 || m1Y==0) 
			bTran0=1;
		if(m1X==cells[1]) 
		{
			m1X=0;
			bTran0=1;
		}
		else if(m1X==cells[1]-1 || m1X==0)
			bTran0=1;
		
		m1Z++;
		m1Y++;
		m1X++;

		//m1=((m1Z-1)*cells[2]+m1Y-1)*cells[1]+m1X+nAtom;

		for (offset=1;offset<=14;offset++)
		{	
			bTran=bTran0;
			m2Z=m1Z+iofZ[offset];
			shift[3]=0.0f;
			if(m2Z>cells[3])
			{
				m2Z=1;
				shift[3]=region[3];bTran=1;				
			}
			else if(m2Z==0)
			{
				m2Z=cells[3];
				shift[3]=-region[3];bTran=1;
			}

			m2X=m1X+iofX[offset];
			shift[1]=0.0f;
			if(m2X>cells[1])
			{
				m2X=1;
				shift[1]=region[1];bTran=1;
			}
			else if(m2X==0)
			{
				m2X=cells[1];
				shift[1]=-region[1];bTran=1;
			}

			m2Y=m1Y+iofY[offset];
			shift[2]=0.0f;
			if(m2Y>cells[2])
			{
				m2Y=1;
				shift[2]=region[2];bTran=1;
			}
			else if(m2Y==0)
			{
				m2Y=cells[2];
				shift[2]=-region[2];bTran=1;
			}
			
			m2=((m2Z-1)*cells[2]+m2Y-1)*cells[1]+m2X+nAtom;
			
			j2=cellList[m2];

			if(offset==1) j2=cellList[j1];

			while(j2>0)
			{
				if(inPolyj1!=0 && inPoly[j2]==inPolyj1)
				{
					j22=j2-ChainBegin[i];
					if(ExcludeMatrix[i][j11%ChainLength[i]][j22%ChainLength[i]]==1)
					goto flag;
				}
				dr[1]=xtmp-rx[j2]-shift[1];
				dr[2]=ytmp-ry[j2]-shift[2];
				dr[3]=ztmp-rz[j2]-shift[3];
				
				rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];
				if(rr<rrNebr && (!(type2[j1] ==1 && type2[j2]==1)) )
				{
					if(bTran==0)
					{
						nebrTabLen++;
#ifdef MYDEBUG
						if(nebrTabLen>=nebrTabMax-2)
							return 1;
#endif
						nebrTab2[nebrTabLen]=j2;
					}
					else
					{
						nebrTabLen2++;
#ifdef MYDEBUG
						if(nebrTabLen2>=nebrTabMax-2)
							return 1;
#endif
						nebrTab4[nebrTabLen2]=j2;
					}
				} // if(rr<rrNebr)
flag:
				j2=cellList[j2];
			}
		}// for(offset...)
	}// for(j1=1;j1<=nAtom;j1++);
	nebrTab1[j1]=nebrTabLen+1;
	nebrTab3[j1]=nebrTabLen2+1;
	return 0;
}


int  BuildNebrList_second() //nebr list method
{
	float dr[4], invWid[4],shift[4];
	float rr,rrNebr;
	int   bTran,bTran0;
	float xtmp,ytmp,ztmp;
	int   c,j1,j2,k,m1X,m1Y,m1Z,m2,m2X,m2Y,m2Z,n,offset;

	// below three array define the half-neighbour cell(total 14), include itself.
	int    iofX[]={0,0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1};
	int    iofY[]={0,0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1};
	int    iofZ[]={0,0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1};
	int    c1,c2,c3;

	nebrNow_second=0;  //clear flag for buildneighborlist
	dispHi_second=0.0f;   //clear the sum of max V, which is used to decide whether to buildlist or not!

	rrNebr=(rCut_second+rNebrShell_second)*(rCut_second+rNebrShell_second);
	for(k=1;k<=3;k++)
		invWid[k]=cells_second[k]/region[k];
//below code: cellList play dual role: The first parts of the array consists of pointers linking
//different atoms belong to the same cell, while the second parts one per cell, point to the
//first atom in each cell, 0 mean the final atom in one list of cell and an empty cell respectively
	for(n=nAtom+1;n<=nAtom+cells_second[1]*cells_second[2]*cells_second[3];n++)
		cellList[n]=0;

	
	for(n=1;n<=nAtom;n++)
	{
		if(type[n]!=5)
			continue;
		c3=(int)((rz[n]+regionH[3])*invWid[3]);
		c2=(int)((ry[n]+regionH[2])*invWid[2]);
		c1=(int)((rx[n]+regionH[1])*invWid[1]);
		if(c3==cells[3]) c3=0;
		if(c2==cells[2]) c2=0;
		if(c1==cells[1]) c1=0;
		c=((int)c3*cells[2]+(int)c2)*cells[1]+(int)c1+nAtom+1;
		if(c>nAtom+cells[1]*cells[2]*cells[3] || c<nAtom+1)
		{
			printf("In step %d for atom number %d \r\n",stepCount,n);
			printf(" 3  2  1 is %d,%d,%d\r\n", cells[3],cells[2],cells[1]);
			printf("c3,c2,c1 is %d,%d,%d\r\n",c3,c2,c1);
			printf("n is %d, rx[n] is %f, regionH[1] is %f, invWid[1] is %f\r\n",n,rx[n],regionH[1],invWid[1]);
			printf("n is %d, ry[n] is %f, regionH[2] is %f, invWid[2] is %f\r\n",n,ry[n],regionH[2],invWid[2]);
			printf("n is %d, rz[n] is %f, regionH[3] is %f, invWid[3] is %f\r\n",n,rz[n],regionH[3],invWid[3]);
			return 2;  //over flow !!!
		}
		cellList[n]=cellList[c];
		cellList[c]=n;
	}
	nebrTabLen_second=0;
	nebrTabLen2_second=0;
	
	for(j1=nstart;j1<=nend;j1++)
	{
		nebrTab1_second[j1]=nebrTabLen_second+1;
		nebrTab3_second[j1]=nebrTabLen2_second+1;

		if(type[j1]!=5)
			continue;


		xtmp=rx[j1];
		ytmp=ry[j1];
		ztmp=rz[j1];

		bTran0=0;
		m1Z=(int)((rz[j1]+regionH[3])*invWid[3]);
		m1Y=(int)((ry[j1]+regionH[2])*invWid[2]);
		m1X=(int)((rx[j1]+regionH[1])*invWid[1]);
		if(m1Z==cells[3]) 
		{
			m1Z=0;
			bTran0=1;
		}
		else if(m1Z==cells[3]-1 || m1Z==0)
			bTran0=1;
		if(m1Y==cells[2])
		{
			m1Y=0;
			bTran0=1;
		}
		else if(m1Y==cells[2]-1 || m1Y==0) 
			bTran0=1;
		if(m1X==cells[1]) 
		{
			m1X=0;
			bTran0=1;
		}
		else if(m1X==cells[1]-1 || m1X==0)
			bTran0=1;
		
		m1Z++;
		m1Y++;
		m1X++;

		//m1=((m1Z-1)*cells[2]+m1Y-1)*cells[1]+m1X+nAtom;

		for (offset=1;offset<=14;offset++)
		{	
			bTran=bTran0;
			m2Z=m1Z+iofZ[offset];
			shift[3]=0.0f;
			if(m2Z>cells[3])
			{
				m2Z=1;
				shift[3]=region[3];bTran=1;				
			}
			else if(m2Z==0)
			{
				m2Z=cells[3];
				shift[3]=-region[3];bTran=1;
			}

			m2X=m1X+iofX[offset];
			shift[1]=0.0f;
			if(m2X>cells[1])
			{
				m2X=1;
				shift[1]=region[1];bTran=1;
			}
			else if(m2X==0)
			{
				m2X=cells[1];
				shift[1]=-region[1];bTran=1;
			}

			m2Y=m1Y+iofY[offset];
			shift[2]=0.0f;
			if(m2Y>cells[2])
			{
				m2Y=1;
				shift[2]=region[2];bTran=1;
			}
			else if(m2Y==0)
			{
				m2Y=cells[2];
				shift[2]=-region[2];bTran=1;
			}
			
			m2=((m2Z-1)*cells[2]+m2Y-1)*cells[1]+m2X+nAtom;
			
			j2=cellList[m2];

			if(offset==1) j2=cellList[j1];

			while(j2>0)
			{
				dr[1]=xtmp-rx[j2]-shift[1];
				dr[2]=ytmp-ry[j2]-shift[2];
				dr[3]=ztmp-rz[j2]-shift[3];
				
				rr=dr[1]*dr[1]+dr[2]*dr[2]+dr[3]*dr[3];
				if(rr<rrNebr)
				{
					if(bTran==0)
					{
						nebrTabLen_second++; 
						nebrTab2_second[nebrTabLen_second]=j2;
					}
					else
					{
						nebrTabLen2_second++; 
						nebrTab4_second[nebrTabLen2_second]=j2;
					}
				} // if(rr<rrNebr)
				j2=cellList[j2];
			}
		}// for(offset...)
	}// for(j1=1;j1<=nAtom;j1++);
	nebrTab1_second[j1]=nebrTabLen_second+1;
	nebrTab3_second[j1]=nebrTabLen2_second+1;
	return 0;
}



void ApplyBoundaryCond()
{
	int n;

	for(n=nstart;n<=nend;n++)
	{
		if(rx[n]>=regionH[1])
		{
			rx[n]-=region[1];
		}
		else if(rx[n]<-regionH[1])
		{
			rx[n]+=region[1];
		}
		
		if(ry[n]>=regionH[2])
		{
			ry[n]-=region[2];
		}
		else if(ry[n]<-regionH[2])
		{
			ry[n]+=region[2];
		}
			
		if(rz[n]>=regionH[3])
		{
			rz[n]-=region[3];
		}
		else if(rz[n]<-regionH[3])
		{
			rz[n]+=region[3];
		}	
	}
}



void EvalProps()
{
}

void  EvalBuildNeighbor()
{
	float vv;
	int n;
	vvMax=0.0f;               
	for(n=1;n<=nAtom;n++)
	{
		vv=rvx[n]*rvx[n]+rvy[n]*rvy[n]+rvz[n]*rvz[n];
		if(vv>vvMax)
			vvMax=vv; 
	}

#ifdef MYMPI
	vv=vvMax;
	MPI_Allreduce(&vv,&vvMax,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
#endif

	dispHi=dispHi+(float)(sqrt(vvMax)*deltaT);    
#ifdef MYDEBUG
	printf("In Step %d, vvMax is %f, and DispHi is %f\r\n",stepCount,vvMax,dispHi);
#endif
	if(dispHi>0.5*rNebrShell) nebrNow=1;
}

void  EvalBuildNeighbor_second()
{
	float vv;
	int n;
	vvMax=0.0f;               
	for(n=1;n<=nAtom;n++)
	{
		if(type[n]!=5)
			continue;
		vv=rvx[n]*rvx[n]+rvy[n]*rvy[n]+rvz[n]*rvz[n];
		if(vv>vvMax)
			vvMax=vv; 
	}

#ifdef MYMPI
	vv=vvMax;
	MPI_Allreduce(&vv,&vvMax,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
#endif

	dispHi_second=dispHi_second+(float)(sqrt(vvMax)*deltaT);    
	if(dispHi_second>0.5*rNebrShell_second) nebrNow_second=1;
}

void AccumProps(int icode)
{
	float vv;
	int k,kk,n;

	if(icode==0)
	{
		kinEnergy=0.0f;
		potEnergy=0.0f;
		hEnergy=0.0f;
		hydroEnergy=0.0f;
		hNumber=0.0f;
		bondEnergy=0.0f;
		angleEnergy=0.0f;
		torsionEnergy=0.0f;

		spressure=0.0f;
		sspressure=0.0f;
		svolume=0.0f;
		ssvolume=0.0f;
		
		for(k=1;k<=3;k++)
		for(kk=1;kk<=3;kk++)
		{
			svirSum[k][kk]=0.0f;
		}
		svirSum[0][0]=0.0f;

	}
	else if(icode==1)
	{
		vvSum=0.0f;
		for(k=1;k<=3;k++)
		for(kk=1;kk<=3;kk++)
		{
			kinSum[k][kk]=0.0f;
		}

		for(n=nstart;n<=nend;n++)
		{
			vv=rvx[n]*rvx[n]+rvy[n]*rvy[n]+rvz[n]*rvz[n];
			vv/=mm[type[n]];
			kinSum[1][1]+=rvx[n]*rvx[n];
			kinSum[1][2]+=rvx[n]*rvy[n];
			kinSum[1][3]+=rvx[n]*rvz[n];
			kinSum[2][1]+=rvy[n]*rvx[n];
			kinSum[2][2]+=rvy[n]*rvy[n];
			kinSum[2][3]+=rvy[n]*rvz[n];
			kinSum[3][1]+=rvz[n]*rvx[n];
			kinSum[3][2]+=rvz[n]*rvy[n];
			kinSum[3][3]+=rvz[n]*rvz[n];
			vvSum=vvSum+vv;
		}  

		vVal=vvSum/nAtom/3;
		
		volume=region[1]*region[2]*region[3];
		pressure=(vvSum+virSum[0][0])/(volume*3);
		
		virSum[0][0]/=(volume*3);
		svirSum[0][0]+=virSum[0][0];

		for(k=1;k<=3;k++)
		for(kk=1;kk<=3;kk++)
		{
			virSum[k][kk]/=volume;
			kinSum[k][kk]/=volume;
			svirSum[k][kk]+=virSum[k][kk];
			skinSum[k][kk]+=kinSum[k][kk];
		}	
		
		if(vVal==0) 
			vVal=0.01f;
		lamad=(float)sqrt(1+deltaT*(coup_temperature/vVal-1));

		spressure+=pressure;
		sspressure+=pressure*pressure;
		volume/=1000;
		svolume+=volume;
		ssvolume+=volume*volume;
		
		kinEnergy+=vVal;
		potEnergy+=uVal/nAtom;
		hEnergy+=hVal/nAtom;
		hNumber+=hnumberVal;
		outofplaneNumber+=outofplaneVal;
		hydroEnergy+=hydroVal/nAtom;
		bondEnergy+=bondVal/nAtom;
		angleEnergy+=angleVal/nAtom;
		torsionEnergy+=torsionVal/nAtom;

	}
	else if(icode==2)
	{
		spressure/=stepAvg;
		sspressure=(float)sqrt(sspressure/stepAvg-spressure*spressure);
		svolume/=stepAvg;
		ssvolume/=stepAvg;
		ssvolume-=svolume*svolume;
		if(ssvolume>=0) 
			ssvolume=(float)sqrt(ssvolume);
		svirSum[0][0]/=stepAvg;

		for(k=1;k<=3;k++)
		for(kk=1;kk<=3;kk++)
		{
			svirSum[k][kk]/=stepAvg;
			skinSum[k][kk]/=stepAvg;
		}

		vvSum=0.0f;			
		for(n=1;n<=nAtom;n++)
		{
			vv=rvx[n]*rvx[n]+rvy[n]*rvy[n]+rvz[n]*rvz[n];
			vv/=mm[type[n]];
			vvSum=vvSum+vv;
		}
		
		kinEnergy       /=stepAvg;		
		potEnergy       /=stepAvg;
		hEnergy         /=stepAvg;
		hydroEnergy     /=stepAvg;
		hNumber         /=stepAvg;
		outofplaneNumber/=stepAvg;
		bondEnergy      /=stepAvg;
		angleEnergy     /=stepAvg;
		torsionEnergy   /=stepAvg;
	}
}

float Integrate(float *f,int nf)
{
	float s;
	int i;
	s=(float)(0.5*(f[1]+f[nf]));
	for(i=2;i<=nf-1;nf++)
		s=s+f[i];
	return s;
}

void AllocateMemory()
{
	int k;


	cell=new int* [4];
	rold=new float* [4];

	 rx=new float[nAtom+1];
	 ry=new float[nAtom+1];
	 rz=new float[nAtom+1];

	 rvx=new float[nAtom+1];
	 rvy=new float[nAtom+1];
	 rvz=new float[nAtom+1];
	 
	 rax=new float[nAtom+1];
	 ray=new float[nAtom+1];
	 raz=new float[nAtom+1];

#ifdef MYMPI
	 raxx=new float[nAtom+1];
	 rayy=new float[nAtom+1];
	 razz=new float[nAtom+1];
#endif

	for(k=1;k<=3;k++)
	{
		cell[k]=new int [nAtom+1];
		rold[k]=new float[nAtom+1];
	}

	    type=new int[nAtom+1];
	   type2=new int[nAtom+1];
	  inPoly=new int[nAtom+1];
	cellList=new int[nAtom+cells[1]*cells[2]*cells[3]]; //for cell method
	nebrTab1=new int[nAtom+2];
	nebrTab2=new int[nebrTabMax];
	nebrTab3=new int[nAtom+2];
	nebrTab4=new int[nebrTabMax];
	
	nebrTab1_second=new int[nAtom+2];
	nebrTab2_second=new int[nebrTabMax];
	nebrTab3_second=new int[nAtom+2];
	nebrTab4_second=new int[nebrTabMax];

#ifdef getdiffusion
	r0=new float *[4];
	cell0=new int *[4];
	for(k=1;k<=3;k++)
	{
		r0[k]=new float[nChain+1];
		cell0[k]=new int[nChain+1];
	}
#endif

#ifdef GETMONOMER
	headneib= new int[nChain+1];
#endif

}

void ReleaseMemory()
{
	int k;

	for(k=1;k<=3;k++)
	{
		delete cell[k];
	}
	delete rx;
	delete ry;
	delete rz;
	delete rvx;
	delete rvy;
	delete rvz;
	delete rax;
	delete ray;
	delete raz;
	delete type;
	delete type2;
	delete inPoly;
	delete cell;

#ifdef MYMPI
	delete raxx;
	delete rayy;
	delete razz;
#endif

	//delete cellList, for nebr list method
	delete nebrTab1;
	delete nebrTab2; 
	delete nebrTab3;
	delete nebrTab4; 

	int ii,i;

	for(ii=1;ii<=TypePoly;ii++)
	{
		for(i=0;i<AtomNumber[ii];i++)
		{
			delete ExcludeMatrix[ii][i];
		}						  
		delete 	ExcludeMatrix[ii];

		delete  BondIndex1[ii];
		delete  BondIndex2[ii];
		delete  BondLength[ii];
		delete  BondStrenth[ii];
		
		delete  BondAngleIndex1[ii];
		delete  BondAngleIndex2[ii];
		delete  BondAngleIndex3[ii];
		delete  BondAngle[ii];
		delete  BondAngleStrenth[ii];

		delete  BondTorsionIndex1[ii];
		delete  BondTorsionIndex2[ii];
		delete  BondTorsionIndex3[ii];
		delete  BondTorsionIndex4[ii];
		delete  BondTorsionType[ii];
		delete  BondTorsionStrenth[ii];

		delete  OutofPlaneIndex1[ii];
		delete  OutofPlaneIndex2[ii];
		delete  OutofPlaneIndex3[ii];
		delete  OutofPlaneIndex4[ii];
		delete  OutofPlaneStrenth[ii];

		delete  Interaction14Index1[ii];
		delete  Interaction14Index2[ii];

		delete  ShakeIndex1[ii];
		delete  ShakeIndex2[ii];
		delete  ShakeLength[ii];
	}
}


//These dininations are used for fuction ran. 
#define IA	16807
#define IM	2147483647
#define AM	(1.0/IM)
#define IQ	127773
#define IR	2836
#define NTAB	32
#define NDIV	(1+(IM-1)/NTAB)
#define EPS		1.2e-7
#define RNMX	(1.-EPS)

//generate random number, peroid is nearly 10 pow 12
double ran(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if( *idum<=0 || !iy)
	{
		if( -(*idum)<1) *idum=1;
		else *idum=-(*idum);
		for(j=NTAB+7; j>=0; j--)
		{
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if(*idum<0) *idum+=IM;
			if(j<NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if(*idum<0) *idum+=IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j]=*idum;
	if((temp=AM*iy)>RNMX) return RNMX;
	else return temp;
}

int ProcessError(int iCode)
{
	switch(iCode)
	{
	case 1:
		if( byteread< nAtom+1)
		{
			printf("Error in reading str file! \r\n");
			fout<<"Error in reading str file! \r\n";
			ReleaseMemory();
			fclose(fp);
			return 0;
		}
		break;
	case 2:
		if(BuildListReturn==1)
		{
			printf("Error! build list overflow when step %d \r\n", stepCount);
			fout<<"Error, build list overflow when "<<stepCount<<"\r\n";
			fout.close();
			return 0;
		}
		else if(BuildListReturn==2)
		{
			printf("Error! Time Step too large when step %d \r\n", stepCount);
			fout<< "Error! Time Step too large when "<<stepCount<<"\r\n";
			fout.close();
			return 0;
		}
		break;
	case 3:
		if(ComputeBondForceReturn==1)
		{
			ReleaseMemory();
			printf("Error, Bond too long when step %d! \r\n",stepCount);
			fout<<"Error, Bond too long when "<<stepCount<<"\r\n";
			fout.close();
			return 0;
		}
		else if(ComputeBondForceReturn==2)
		{
			ReleaseMemory();
			printf("Error, Bond too small when step %d \r\n",stepCount);
			fout<<"Error, Bond too small when "<<stepCount<<"\r\n";
			fout.close();
			return 0;
		}
		break;
	case 4:
		if(LeapfrogStepReturn==1) 
		{	
			OutPut(0);
			printf("Error in update r,rv \r\n");
			fout<<"Error in update r,rv \r\n";
			fout.close();
			ReleaseMemory();
			return 0;
		}
		break;
	}
	return 1;
}


void Outputdata()
{
	sprintf(use,"%3d, %5.3f, %5.3f, ",(stepCount+beginnumber)/10000,kinEnergy,potEnergy);
	myprintf(use);
	sprintf(use,"%6.4f, %6.2f, %6.4f, %6.3f, %6.3f, %6.3f, %6.4f",hEnergy,hNumber,hydroEnergy,bondEnergy,angleEnergy,torsionEnergy,outofplaneNumber);
	myprintf(use);
	if(node==0) printf("\r\n");
	sprintf(use,", %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f",
		svirSum[1][1],skinSum[1][1],svirSum[2][2],skinSum[2][2],svirSum[3][3],skinSum[3][3]);
	fout<<use;
	fout<<"\r\n";
	fout.flush();

	if(node==0) OutPut(stepCount+beginnumber);

}

void OutPut(int Index)
{
	char string[255];
	FILE *fp;

	if(Index==0)
		sprintf(string,"%sdpddamp.str",maindirectory);
	else if(Index==-1)
		sprintf(string,"%smin.str",maindirectory);
	else
		sprintf(string,"%ss%d.str",maindirectory,Index);
	fp=fopen(string,"w+b");

	fwrite(&density,sizeof(float),1,fp);
	fwrite(&TypePoly,sizeof(int),1,fp);
	fwrite(nChain+1,sizeof(int),TypePoly,fp);
	fwrite(ChainLength+1,sizeof(int),TypePoly,fp);
	fwrite(ChainBegin+1,sizeof(int),TypePoly+1,fp);
	fwrite(&nAtom,sizeof(int),1,fp);
	fwrite(region+1,sizeof(float),3,fp);
	
	fwrite(type,sizeof(int),(nAtom+1),fp);
	fwrite(type2,sizeof(int),(nAtom+1),fp);
	fwrite(inPoly,sizeof(int),(nAtom+1),fp);
	fwrite(cell[1],sizeof(int),nAtom+1,fp);
	fwrite(cell[2],sizeof(int),nAtom+1,fp);
	fwrite(cell[3],sizeof(int),nAtom+1,fp);
	fwrite(rx,sizeof(float),(nAtom+1),fp);
	fwrite(ry,sizeof(float),(nAtom+1),fp);
	fwrite(rz,sizeof(float),(nAtom+1),fp);
	fwrite(rvx,sizeof(float),(nAtom+1),fp);
	fwrite(rvy,sizeof(float),(nAtom+1),fp);
	fwrite(rvz,sizeof(float),(nAtom+1),fp);
	fclose(fp);
}

void Receivecommand()
{
	ifstream fmes(messfilename);
	fmes>>use;
	if(strcmp(use,"end")==0)
	{
		moreCycles=0;
	}
	fmes.close();
}

void myprintf(char *ps)
{
	if(node==0) printf(ps);
	fout<<ps;
}


int Setparameter()
{
	int k;
	int i,j,ii;

	//open file for output data and message!
	sprintf(use,"data%d.txt",node);
	fout.open(use);  

	//read force parameter from init file.
	sprintf(use,"read from init file %s.\r\n",initfilename);
	myprintf(use);

	ifstream fin(initfilename);
	fin>>use>>version;
	if( (version/10)<30 )
	{
		myprintf("Version of init don't(dosn't) coincide\r\n");
		return 0;
	}
	fin>>use>>minimizer;
	fin>>use>>Min_Step;         //basic dissp force
	fin>>use>>deltaT;
	fin>>use>>stepAvg;
	fin>>use>>stepLimit;       
	fin>>use>>nebrTabFac;      // length of neighbour list
	fin>>use>>rNebrShell;      // shell= this value+rcut
	fin>>use>>coup_temperature;          

	fin>>use>>coup_pressure;
	fin>>use>>Langevin;
	fin>>use>>Langevin;
	fin>>use>>Langevin;
	fin>>use>>Langevin;
	fin>>use>>Langevin;
	
	fin>>use>>constraintPrec;


	fin>>use>>beginnumber;
	fin>>use>>seed;
	fin>>use>>seed;
	fin>>use>>beginnumber;
	fin.close();

	if(beginnumber<0)
	{
		myprintf("Something is wrong in ini file! \r\n");
		return 0;
	}

	sprintf(use,"\r\ndelta T    = %f\r\n",deltaT);	
	myprintf(use);
#ifdef SSL
	sprintf(use,"!!!!!! Using SSE !!!!!!!\r\n");
	myprintf(use);
#endif

#ifndef BONDANGLE
	sprintf(use,"Note that Bond Angle force is NOT included!\r\n");
	myprintf(use);
#endif

#ifdef PRESSURECOUP
#ifdef ANISOTROPIC
	sprintf(use,"Applying anisotropic pressure coup! pressure = %f\r\n",coup_pressure);
#else
	sprintf(use,"Applying isotropic pressure coup! pressure = %f\r\n",coup_pressure);
#endif
#else 
	sprintf(use,"Constant volume simulation ! \r\n");
#endif
	myprintf(use);

#ifdef TEMPERATURECOUP
	sprintf(use,"Applying temeprature coup ! temp = %f\r\n",coup_temperature);
	myprintf(use);
#endif

#ifdef LANGEVIN
	sprintf(use,"Applying langevin dynamics ! temp = %f\r\n",coup_temperature);
	myprintf(use);
#endif

#ifndef TEMPERATURECOUP
#ifndef LANGEVIN
	sprintf(use,"Constant energy simulation ! \r\n");
	myprintf(use);
#endif
#endif

	sprintf(use,"Total %d Steps, output every %d steps\r\n",stepLimit,stepAvg);
	myprintf(use);

	
	//*******************************************************************************
	//read cell parameter from structures file.
	sprintf(use,"%ss%d.str",maindirectory,beginnumber);
	fp=fopen(use,"r+b");

	if(fp==NULL)
	{
		printf("structure file %s is not found! \r\n",use);
		fout <<"structure  file "<<use<<" is not found! \r\n";
		return 0;
	}
	fread(&density,sizeof(float),1,fp);
	fread(&TypePoly,sizeof(int),1,fp);
	fread(nChain+1,sizeof(int),TypePoly,fp);
	fread(ChainLength+1,sizeof(int),TypePoly,fp);
	fread(ChainBegin+1,sizeof(int),TypePoly+1,fp);
	fread(&nAtom,sizeof(int),1,fp);
	fread(region+1,sizeof(float),3,fp);
	

	sprintf(use,"read from strucures file: %ss%d.str\r\n",maindirectory,beginnumber);	
	myprintf(use);

#ifdef MYMPI
	nstart= (int)(((float)node)/nprocs*nAtom + 1);
    nend  = (int)(((float)node+1)/nprocs*nAtom);
    nlocal= nend - nstart + 1;
#else
	nstart=1;
	nend=nAtom;
#endif

	sprintf(use,"I am %d, handle from %d to %d of total %d\n\r",node,nstart,nend,nAtom);
	myprintf(use);

	sprintf(use,"cell edge is: %f * %f * %f, volume is %f\r\n",region[1],region[2],region[3],region[1]*region[2]*region[3]);	
	myprintf(use);

#ifdef MYMPI
	for(i=1;i<=TypePoly;i++)
	{
		nchainstart[i] = (int)(((float)node)/nprocs*nChain[i] + 1);
		nchainend[i]   = (int)(((float)node+1)/nprocs*nChain[i]);
	}

	for(i=1;i<=TypePoly;i++)
	{
		sprintf(use,"Chain length is %d, Number is %d\r\n",ChainLength[i],nChain[i]);
		myprintf(use);
		sprintf(use, "I am %d, handle from %d to %d\n\r",node,nchainstart[i],nchainend[i]);
		myprintf(use);
	}
#else
	for(i=1;i<=TypePoly;i++)
	{
		sprintf(use,"Chain length is %d, Number is %d\r\n",ChainLength[i],nChain[i]);
		myprintf(use);
	}
#endif

	//Set some params
	ddeltaTH=deltaT*deltaT/2;
	sqrtdeltaT=(float)(sqrt(deltaT));
	pi=(float)(4.0*atan(1.0));

	rCut=4.5f/2.0f;          
	rrCut=rCut*rCut; 
	rCut_second=8.0f/2.0f;
	rrCut_second=rCut_second*rCut_second;

	rNebrShell_second=rNebrShell*2; //?? should tune for the best value !!
	
	nebrTabMax=nAtom*nebrTabFac;					 //neighbor list
	
	for(k=1;k<=3;k++)							 // this loop for cells method only
	{
		cells[k]=(int)(region[k]/(rCut+rNebrShell)); //+rNebrShell used for nerghbor list method
		cells_second[k]=(int)(region[k]/(rCut_second+rNebrShell_second)); //+rNebrShell used for nerghbor list method
		regionH[k]=region[k]/2;
	}
	
	dim=3;
	
	//read structures
	AllocateMemory();
	fread(type,sizeof(int),(nAtom+1),fp);
	fread(type2,sizeof(int),(nAtom+1),fp);
	fread(inPoly,sizeof(int),(nAtom+1),fp);
	fread(cell[1],sizeof(int),nAtom+1,fp);
	fread(cell[2],sizeof(int),nAtom+1,fp);
	fread(cell[3],sizeof(int),nAtom+1,fp);
	fread(rx,sizeof(float),(nAtom+1),fp);
	fread(ry,sizeof(float),(nAtom+1),fp);
	fread(rz,sizeof(float),(nAtom+1),fp);
	fread(rvx,sizeof(float),(nAtom+1),fp);
	fread(rvy,sizeof(float),(nAtom+1),fp);
	byteread=fread(rvz,sizeof(float),(nAtom+1),fp);
	if(ProcessError(1)==0) return 0;
	fclose(fp);



//read from force field and topological file
	int m,n;
	float temp_delta;
	ifstream fin2("test.top");

	fin2>>TypeAtom;
	if(TypeAtom>=TYPE)
	{
		sprintf(use,"Error ! macro TYPE is too small !, should be at least %d",TypeAtom+1);
		myprintf(use);
		return 0;
	}
	for(i=1;i<=TypeAtom;i++)
	{
		fin2>>mm[i];
		mm[i]=1.0f/mm[i];
	}

	fin2>>ForceNumber;
	for(i=1;i<=ForceNumber;i++)
	{
		fin2>>m;
		fin2>>n;
		fin2>>vmForceMatrix[m][n];
		fin2>>vmForceMatrix14[m][n]; 
		fin2>>temp_delta;
		temp_delta*=temp_delta;
		temp_delta*=temp_delta;
		temp_delta=temp_delta*temp_delta*temp_delta;
		vmDeltaMatrix[m][n]=temp_delta;
		fin2>>temp_delta;
		temp_delta*=temp_delta;
		temp_delta*=temp_delta;
		temp_delta=temp_delta*temp_delta*temp_delta;
		vmDeltaMatrix14[m][n]=temp_delta;

		vmForceMatrix  [n][m]=vmForceMatrix  [m][n];
		vmForceMatrix14[n][m]=vmForceMatrix14[m][n];
		vmDeltaMatrix  [n][m]=vmDeltaMatrix  [m][n];
		vmDeltaMatrix14[n][m]=vmDeltaMatrix14[m][n];
	}
	fin2>>H_Force;
	fin2>>Hydro_Force; 

	fin2>>TypePoly;
	if(TypePoly>=TYPEPOLY)
	{
		sprintf(use,"Error ! macro TYPEPOLY is too small !, should be at least %d",TypePoly+1);
		myprintf(use);
		return 0;
	}

	sprintf(use,"Total %d types of Atoms, %d types of Protein.\n",TypeAtom,TypePoly);
	myprintf(use);
	sprintf(use,"Hydrogen force is %f and hydrophobic force is %f\n",H_Force,Hydro_Force);
	myprintf(use);

	for(ii=1;ii<=TypePoly;ii++)
	{

		fin2>>AtomNumber[ii];
		ExcludeMatrix[ii]=new int * [AtomNumber[ii]];
		for(i=0;i<AtomNumber[ii];i++)
		{
			ExcludeMatrix[ii][i]=new int[AtomNumber[ii]];
		}		
		for(i=0;i<AtomNumber[ii];i++)
		for(j=0;j<AtomNumber[ii];j++)
			ExcludeMatrix[ii][i][j]=0;
		printf("succeed in exclude...");
		
		fin2>>BondNumber[ii];
		BondIndex1[ii]=new int[BondNumber[ii]+1];
		BondIndex2[ii]=new int[BondNumber[ii]+1];
		BondLength[ii]=new float[BondNumber[ii]+1];
		BondStrenth[ii]=new float[BondNumber[ii]+1];
		for(i=1;i<=BondNumber[ii];i++)
		{
			fin2>>BondIndex1[ii][i];
			fin2>>BondIndex2[ii][i];
			fin2>>BondLength[ii][i];
			fin2>>BondStrenth[ii][i];
			m=BondIndex1[ii][i];
			n=BondIndex2[ii][i];
			if(m==ChainLength[ii])
				m=0;
			if(n==ChainLength[ii])
				n=0;
			ExcludeMatrix[ii][m][n]=1;
			ExcludeMatrix[ii][n][m]=1;
		}	
		printf("Bond...");

		
		fin2>>BondAngleNumber[ii];
		BondAngleIndex1[ii]=new int[BondAngleNumber[ii]+1];
		BondAngleIndex2[ii]=new int[BondAngleNumber[ii]+1];
		BondAngleIndex3[ii]=new int[BondAngleNumber[ii]+1];
		BondAngle[ii]=new float[BondAngleNumber[ii]+1];
		BondAngleStrenth[ii]=new float[BondAngleNumber[ii]+1];
		for(i=1;i<=BondAngleNumber[ii];i++)
		{
			fin2>>BondAngleIndex1[ii][i];
			fin2>>BondAngleIndex2[ii][i];
			fin2>>BondAngleIndex3[ii][i];
			fin2>>BondAngle[ii][i];
			fin2>>BondAngleStrenth[ii][i];
			BondAngle[ii][i]=BondAngle[ii][i]*pi/180.0f; 
			m=BondAngleIndex1[ii][i];
			n=BondAngleIndex3[ii][i];
			if(m==ChainLength[ii])
				m=0;
			if(n==ChainLength[ii])
				n=0;
			ExcludeMatrix[ii][m][n]=1;
			ExcludeMatrix[ii][n][m]=1;
		}
		printf("Angle...");

		
		fin2>>BondTorsionNumber[ii];
		BondTorsionIndex1[ii]=new int[BondAngleNumber[ii]+1];
		BondTorsionIndex2[ii]=new int[BondAngleNumber[ii]+1];
		BondTorsionIndex3[ii]=new int[BondAngleNumber[ii]+1];
		BondTorsionIndex4[ii]=new int[BondAngleNumber[ii]+1];
		BondTorsionType[ii]=new int[BondAngleNumber[ii]+1];
		BondTorsionStrenth[ii]=new float[BondAngleNumber[ii]+1];
		for(i=1;i<=BondTorsionNumber[ii];i++)
		{
			fin2>>BondTorsionIndex1[ii][i];
			fin2>>BondTorsionIndex2[ii][i];
			fin2>>BondTorsionIndex3[ii][i];
			fin2>>BondTorsionIndex4[ii][i];
			fin2>>BondTorsionType[ii][i];
			fin2>>BondTorsionStrenth[ii][i];
		}
		printf("Torsion...");

		
		fin2>>OutofPlaneNumber[ii];
		OutofPlaneIndex1[ii]=new int[OutofPlaneNumber[ii]+1];
		OutofPlaneIndex2[ii]=new int[OutofPlaneNumber[ii]+1];
		OutofPlaneIndex3[ii]=new int[OutofPlaneNumber[ii]+1];
		OutofPlaneIndex4[ii]=new int[OutofPlaneNumber[ii]+1];
		OutofPlaneStrenth[ii]=new float[OutofPlaneNumber[ii]+1];
		for(i=1;i<=OutofPlaneNumber[ii];i++)
		{
			fin2>>OutofPlaneIndex1[ii][i];
			fin2>>OutofPlaneIndex2[ii][i];
			fin2>>OutofPlaneIndex3[ii][i];
			fin2>>OutofPlaneIndex4[ii][i];
			fin2>>OutofPlaneStrenth[ii][i];
		}
 		printf("Out_P...");

		fin2>>Interaction14Number[ii];
		Interaction14Index1[ii]=new int[Interaction14Number[ii]+1];
		Interaction14Index2[ii]=new int[Interaction14Number[ii]+1];  
		for(i=1;i<=Interaction14Number[ii];i++)
		{
			fin2>>Interaction14Index1[ii][i];
			fin2>>Interaction14Index2[ii][i];
			m=Interaction14Index1[ii][i];
			n=Interaction14Index2[ii][i];
			if(m==ChainLength[ii])
				m=0;
			if(n==ChainLength[ii])
				n=0; 
			ExcludeMatrix[ii][m][n]=1;
			ExcludeMatrix[ii][n][m]=1;
		}
		printf("14...");


		fin2>>ShakeNumber[ii];
		ShakeIndex1[ii]=new int[ShakeNumber[ii]+1];
		ShakeIndex2[ii]=new int[ShakeNumber[ii]+1];
		ShakeLength[ii]=new float[ShakeNumber[ii]+1];
		for(i=1;i<=ShakeNumber[ii];i++)
		{
			fin2>>ShakeIndex1[ii][i];
			fin2>>ShakeIndex2[ii][i];
			fin2>>ShakeLength[ii][i];
			ShakeLength[ii][i]*=ShakeLength[ii][i];
			
			m=ShakeIndex1[ii][i];
			n=ShakeIndex2[ii][i];
			if(m==ChainLength[ii])
				m=0;
			if(n==ChainLength[ii])
				n=0;
			ExcludeMatrix[ii][m][n]=1;
			ExcludeMatrix[ii][n][m]=1;
		}
	}
	printf("Shake...\n");
	
	fin2.close();		 
	
	moreCycles=1;
	return 1;
}