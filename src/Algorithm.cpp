#include"Algorithm.h"
#include<cmath>
#include"VecMat.h"

Lsq::Lsq(const Matrix &Pxx_,const Matrix &L,const Matrix &H,int n)
{
    l=L;
    Matrix mixed(Pxx_.get_row()+H.get_row(),n);
    // std::cout<<Pxx_<<'\n'<<Inv_LU(Pxx_)<<'\n'<<Inv_QR(Pxx_)<<std::endl;
    Matrix U(Trans(Cholesky(Inv_LU(Pxx_))));
    // std::cout<<U<<std::endl;
    for(int i=0;i<mixed.get_row();i++)
    {
        for(int j=0;j<n;j++)
        {
            if(i<n)
            {
                mixed(i,j)=U(i,j);
            }
            else
            {
                mixed(i,j)=H(i-n,j);
            }
        }
    }
    B=H;
    // std::cout<<mixed;
    Matrix *M;
    // std::cout<<mixed<<std::endl;
    M=QR_decompose(mixed);
    Matrix temp1(mixed.get_row(),1);
    for(int k=0;k<mixed.get_row();k++)
    {
        if(k<n)
        {
            temp1(k,0)=0.0;
        }
        else
        {
            temp1(k,0)=l(k-n,0);
        }
    }
    // std::cout<<l<<M[0]<<std::endl;
    Matrix temp(Trans(M[0])*temp1);
    Matrix X(n,1);
    double sum;
    X(n-1,0)=temp(n-1,0)/M[1](n-1,n-1);
    for(int i=n-2;i>=0;i--)
    {
        sum=0.0;
        for(int j=i+1;j<n;j++)
        {
            sum+=M[1](i,j)*X(j,0);
        }
        X(i,0)=(temp(i,0)-sum)/M[1](i,i);
    }
    x=X;
    // std::cout<<x-Inv_LU(Trans(mixed)*mixed)*Trans(mixed)*temp1<<std::endl
    // std::cout<<x<<std::endl;
    v=B*x-l;
    // std::cout<<v<<std::endl;
    
    if(l.get_row()==n)
    {
        std::cerr<<"ERROR:The number of Observations is n,there is no corr !!! "<<std::endl;
    }
    corr=sqrt((Trans(v)*v)(0,0)/(l.get_row()-n));
    
    // Pxx_post=Inv_LU(Trans(B)*B);
    // std::cout<<M[1].slice(0,n-1,0,n-1)<<std::endl;
    Pxx_post=LUinv(M[1].slice(0,n-1,0,n-1))*Trans(LUinv(M[1].slice(0,n-1,0,n-1)));
    
    delete [] M;
}


Lsq::Lsq(const Matrix &L,const Matrix &H,int n)
{
    l=L;
    B=H;
    Matrix *M;
    M=QR_decompose(B);
    Matrix temp(Trans(M[0])*l);
    Matrix X(n,1);
    double sum;
    X(n-1,0)=temp(n-1,0)/M[1](n-1,n-1);
    for(int i=n-2;i>=0;i--)
    {
        sum=0.0;
        for(int j=i+1;j<n;j++)
        {
            sum+=M[1](i,j)*X(j,0);
        }
        X(i,0)=(temp(i,0)-sum)/M[1](i,i);
    }
    x=X;
    v=B*x-l;
    if(l.get_row()==n)
    {
        std::cerr<<"ERROR:The number of Observations is n,there is no corr !!! "<<std::endl;
    } 
    corr=sqrt((Trans(v)*v)(0,0)/(l.get_row()-n));
    // std::cout<<M[1]<<std::endl;
    Pxx_post=LUinv(M[1].slice(0,n-1,0,n-1))*Trans(LUinv(M[1].slice(0,n-1,0,n-1)));
    // std::cout<<LUinv(M[1].slice(0,n-1,0,n-1))<<'\n'<<Pxx_post<<std::endl;
    delete [] M;
}


void Split(std::string str,const char sign,double *data)
{
    int len,count=0,brepot=0;
    int star_end[2]={0};
    len=str.length();
    for(int i=1;i<len;i++)
    {
        if((str[i]==sign&&str[i-1]!=sign)||i==len-1)
        {
            brepot++;
            star_end[0]=star_end[1];
            star_end[1]=i;
            if(star_end[1]>star_end[0]&&brepot>1)
            {   
                data[count]=atof(str.substr(star_end[0],star_end[1]).c_str());
                star_end[0]=star_end[1];
                count++;
            }
        }
    }
}

