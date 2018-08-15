#include<stdio.h>
#include<string.h>
#include<math.h>
#include<algorithm>
#include<Rinternals.h>
#include <Rdefines.h>
#include<R.h>
//#include"myMatrix2tmpCppCore.h"
#define getDims(A) INTEGER(coerceVector(getAttrib(A,R_DimSymbol),INTSXP))
using namespace std;

extern"C"{
	

void deal(double*,double*sum,double*rsum,int);
void GetSum(double& SumStmpIij,double& SumIij,double* stmp,double *sum,double *rsum,int LenStmp,double First,double Second);
void display(double* X,int len,char* str);
void GetMean(double* X,double* mean,int n,int p);
void GetStd(double* X,double* std,int n,int p);

int n,p,c;


void From8to9(double* X)
{
    double *mean=(double*)malloc(p*sizeof(double));
    GetMean(X,mean,n,p);
    for(int i=0;i<p;++i)
    {
        for(int j=0;j<n;++j)
        {
            *(X+j+i*n)-=mean[i];
        }
    }

    double coef=sqrt((double)n/(n-1));
    double *std;
    std=(double*)malloc(p*sizeof(double));
    GetStd(X,std,n,p);
    for(int i=0;i<p;++i)
    {
        for(int j=0;j<n;++j)
        {
            *(X+j+i*n)=*(X+j+i*n)/std[i]*coef;
        }
    }
    free(mean);
    free(std);
}

void From11to12(double* sd,double* S0)
{
	for(int i=0;i<p;++i){
        for(int j=0;j<p;++j){
            if(i==j){
                *(sd+j+i*p)=1.0;
                *(S0+j+i*p)=-99.0;
            }
            else{
                *(sd+j+i*p)=0;
                *(S0+j+i*p)=-100.0;
            }
        }
    }
}

void From13to22(double* X,double* SS,double* S0,double* sd)
{  
    double* e=(double*)malloc(sizeof(double)*n),mij,sum,sum2;
    for(int j=1;j<p;++j)
    {
        for(int i=0;i<j;++i)
        {
            sum=0;
            for(int k=0;k<n;++k)
            {
                e[k]=*(X+k+i*n)*(*(X+k+j*n));
                sum+=e[k];
            }
            mij=sum/n;
            *(S0+i+j*p)=mij;
            sum2=0;
            for(int k=0;k<n;++k)
            {
                sum2+=(e[k]-mij)*(e[k]-mij);
            }
            *(sd+i+j*p)=sqrt(sum2/(n-1));
        }
    }
    double coef=sqrt((double)n);
    for(int i=0;i<p;++i)
        for(int j=0;j<p;++j)
           *(sd+i*p+j)/=coef;
    free(e);
}



void From24to37(double* X,double* C,double* SS,double* S0,double* sd)
{

    double* Csd =(double*)malloc(p*p*sizeof(double));
    double* stmp=(double*)malloc((p*p+4)*sizeof(double));
    int LenStmp=0;
    
    for(int i=0;i<p*p;++i)
      if(*(S0+i)>-10)
        stmp[LenStmp++]=*(S0+i);
    for(int i=0;i<p*p*c;++i)
        *(SS+i)=1;
    
    double*  sum=(double*)malloc(LenStmp*sizeof(double));
    double* rsum=(double*)malloc(LenStmp*sizeof(double));
    deal(stmp,sum,rsum,LenStmp);
    
    for(int ic=0;ic<c;++ic)
    {
        for(int i=0;i<p;++i)
            for(int j=0;j<p;++j)
             *(Csd+j+i*p)=*(sd+j+i*p)*C[ic];
        
        double First,Second;
        int FirstPos,SecondPos;
        for(int j=1;j<p;++j)
        {
            for(int i=0;i<j;++i)
            {
                double SumStmpIij=0,SumIij=0; 
                if(*(Csd+i+j*p)<-0.0)
                   SumIij=SumStmpIij=0.0;
                else 
                {
                    First=-*(Csd+i+j*p)+*(S0+i+j*p);
                    Second=*(Csd+i+j*p)+*(S0+i+j*p);
                    GetSum(SumStmpIij,SumIij,stmp,sum,rsum,LenStmp,First,Second);                   
                }
                *(SS+i+j*p+ic*p*p)=*(SS+j+i*p+ic*p*p)=SumStmpIij/SumIij;
            }
        }
    }
    free(Csd);free(stmp);free(sum);free(rsum);
}


void AC(double* X,double* C,double* SS,double* S0,int* npc)
{
	n=npc[0];p=npc[1];c=npc[2];
    double* sd=(double*)malloc((p*p+4)*sizeof(double));
	From8to9(X);
	From11to12(sd,S0);
    From13to22(X,SS,S0,sd);
    From24to37(X,C,SS,S0,sd);
    free(sd);
}




	
SEXP myMatrix2tmpC(SEXP X, SEXP C)
{
	int *dimX,*dimC;
	
	double* xptr,*Cptr;
	dimX=getDims(X);
	dimC=getDims(C);
	PROTECT(X=coerceVector(X,REALSXP));
	xptr=REAL(X);
	PROTECT(C=coerceVector(C,REALSXP));
	Cptr=REAL(C);
	
	int npc[3],n,p,c;
	n=dimX[0];p=dimX[1];
	c=dimC[0]*dimC[1];
	npc[0]=n;npc[1]=p;npc[2]=c;
	
	double *pS0,*pSS;
	SEXP S0,SS;	
	PROTECT(S0=NEW_NUMERIC(p*p));
	PROTECT(SS=NEW_NUMERIC(p*p*c));
	pS0=NUMERIC_POINTER(S0);pSS=NUMERIC_POINTER(SS);
	
	AC(xptr,Cptr,pSS,pS0,npc);
	
	SEXP ans,namelist;
	char *names[2]={"SS","S0"};
	PROTECT(namelist=allocVector(STRSXP,2));
	for(int i=0;i<2;++i)
	  SET_STRING_ELT(namelist,i,mkChar(names[i])); 
	
	PROTECT(ans=allocVector(VECSXP,2));
	SET_VECTOR_ELT(ans,0,SS);
	SET_VECTOR_ELT(ans,1,S0);
	setAttrib(ans, R_NamesSymbol, namelist); 
	
	UNPROTECT(6);

	return(ans);

}






void deal(double* stmp,double* sum,double* rsum,int LenStmp)
{
    sort(stmp,stmp+LenStmp);
    rsum[LenStmp-1]=stmp[LenStmp-1];
    sum[0]=stmp[0];
    for(int i=1;i<LenStmp;++i)
    {
        sum[i]=sum[i-1]+stmp[i];
        rsum[LenStmp-i-1]=rsum[LenStmp-i]+stmp[LenStmp-i-1];
    }
}

void GetSum(double& SumStmpIij,double& SumIij,double* stmp,
double *sum,double *rsum,int LenStmp,double First,double Second)
{
    int FirstPos,SecondPos;
    int pl=-1,pr=LenStmp-1;
    int mid;
    while(pl<pr)
    {
        mid=(pl+pr+1)/2;
        if(stmp[mid]<Second) pl=mid;
        else pr=mid-1;
    }
    SecondPos=pl;
    pl=0;pr=LenStmp;
    while(pl<pr)
    {
        mid=(pl+pr)/2;
        if(stmp[mid]>First) pr=mid;
        else pl=mid+1;
    }
    FirstPos=pl;
    if(FirstPos!=LenStmp&&SecondPos!=-1&&FirstPos<=SecondPos)
    {
        SumIij=SecondPos-FirstPos+1;
        SumStmpIij=rsum[FirstPos]+sum[SecondPos]-rsum[0];
    }
    else SumIij=SumStmpIij=0.0;
}

void GetMean(double* X,double* mean,int n,int p)
{
    for(int i=0;i<p;++i)
    {
        double sum=0;
        for(int j=0;j<n;++j)
        {
            sum+=*(X+j+i*n);
        }
        mean[i]=sum/n;
    }
}

void GetStd(double* X,double* std,int n,int p)
{
    double* mean=(double*)malloc(p*sizeof(double));
    GetMean(X,mean,n,p);
    for(int i=0;i<p;++i)
    {
        double sum2=0;
        for(int j=0;j<n;++j)
        {
            sum2+=(*(X+j+i*n)-mean[i])*(*(X+j+i*n)-mean[i]);
        }
        std[i]=sqrt(sum2/(n-1));
    }
    free(mean);
}


}