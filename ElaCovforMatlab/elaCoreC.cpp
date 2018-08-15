#include <yvals.h>
#if (_MSC_VER >= 1600)
#define __STDC_UTF_16__
#endif 
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<algorithm>
#include "mex.h"
using namespace std;

int n,p,c;

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
        std[i]=sqrt((double)sum2/(n-1));
    }
}

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
}

void From11to22(double* X,double* SS,double* S0,double* sd)
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
            *(sd+i+j*p)=sqrt((double)sum2/(n-1));
        }
    }
    double coef=sqrt((double)n);
    for(int i=0;i<p;++i)
        for(int j=0;j<p;++j)
           *(sd+i*p+j)/=coef;
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

void GetSum(double& SumStmpIij,double& SumIij,double* stmp,double *sum,double *rsum,int LenStmp,double First,double Second)
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
}



void AC(double* X,double* C,double* SS,double* S0)
{
    double* sd=(double*)malloc(p*p*sizeof(double));
    From8to9(X);
    From11to22(X,SS,S0,sd);
    From24to37(X,C,SS,S0,sd);
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    if(nrhs!=2) mexErrMsgTxt("Two inputs required.");

    double *X,*C,*SS,*S0;
    int c1,c2;
    n=mxGetM(prhs[0]);
    p=mxGetN(prhs[0]);
    c1=mxGetN(prhs[1]);c2=mxGetM(prhs[1]);
    c=c1>c2?c1:c2;

    X=mxGetPr(prhs[0]);
    C=mxGetPr(prhs[1]);

    plhs[0]=mxCreateDoubleMatrix(1,p*p*c,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(1,p*p,mxREAL);
    SS=mxGetPr(plhs[0]);
    S0=mxGetPr(plhs[1]);

    AC(X,C,SS,S0);
}
