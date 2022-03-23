#include"../include/Algorithm.h"
#include<cmath>
#include"../include/Const.h"
#include"../include/VecMat.h"

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
    Matrix Q,R;
    // std::cout<<mixed<<std::endl;
    QR_decompose(mixed,Q,R);
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
    Matrix temp(Trans(Q)*temp1);
    Matrix X(n,1);
    double sum;
    X(n-1,0)=temp(n-1,0)/R(n-1,n-1);
    for(int i=n-2;i>=0;i--)
    {
        sum=0.0;
        for(int j=i+1;j<n;j++)
        {
            sum+=R(i,j)*X(j,0);
        }
        X(i,0)=(temp(i,0)-sum)/R(i,i);
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
    Pxx_post=LUinv(R.slice(0,n-1,0,n-1))*Trans(LUinv(R.slice(0,n-1,0,n-1)));
    
    // delete [] M;
}


Lsq::Lsq(const Matrix &L,const Matrix &H,int n)
{
    l=L;
    B=H;
    Matrix Q,R;
    QR_decompose(B,Q,R);
    Matrix temp(Trans(Q)*l);
    Matrix X(n,1);
    double sum;
    X(n-1,0)=temp(n-1,0)/R(n-1,n-1);
    for(int i=n-2;i>=0;i--)
    {
        sum=0.0;
        for(int j=i+1;j<n;j++)
        {
            sum+=R(i,j)*X(j,0);
        }
        X(i,0)=(temp(i,0)-sum)/R(i,i);
    }
    x=X;
    v=B*x-l;
    if(l.get_row()==n)
    {
        std::cerr<<"ERROR:The number of Observations is n,there is no corr !!! "<<std::endl;
    } 
    corr=sqrt((Trans(v)*v)(0,0)/(l.get_row()-n));
    // std::cout<<M[1]<<std::endl;
    Pxx_post=LUinv(R.slice(0,n-1,0,n-1))*Trans(LUinv(R.slice(0,n-1,0,n-1)));
    // std::cout<<LUinv(M[1].slice(0,n-1,0,n-1))<<'\n'<<Pxx_post<<std::endl;
    // delete [] M;
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

bool gross_error_detection(int freedom,double sigma,const Matrix &v,const Matrix &Weight)
{
    if(freedom<15)
    {
        std::cerr<<"ERROR:Without T that is smaller than 15"<<std::endl;
        exit(1);
    }

    if((Trans(v)*Weight*v/(sigma*sigma))(0,0)<chi2inv[freedom-15])
    {
        return false;
    }
    else
    {
        return true;
    }
}

void get_coef_mat(const Matrix &obs,const Matrix &ref,double clk,Matrix &H,Matrix &l)
{
    int num=obs.get_row();
    Matrix AuxH(num,4),AuxL(num,1);
    double S0;
    for(int i=0;i<num;i++)
    {
        S0=sqrt(pow(ref(0,0)-obs(i,0),2)+pow(ref(1,0)-obs(i,1),2)+pow(ref(2,0)-obs(i,2),2))+clk;
        AuxL(i,0)=obs(i,3)-S0;
        AuxH(i,0)=(ref(0,0)-obs(i,0))/S0;AuxH(i,1)=(ref(1,0)-obs(i,1))/S0;AuxH(i,2)=(ref(2,0)-obs(i,2))/S0;AuxH(i,3)=1.0;
    }
    H=AuxH;
    l=AuxL;
}

EKF::EKF(int num):epoch(1),time(0.0),n(num){}

void EKF::init(const Matrix &x,const Matrix &P)
{
    X=x;
    Pxx=P;
}

void EKF::init(int epoch_,double time_ ,const Matrix &x,const Matrix &P)
{
    epoch=epoch_;
    time=time_;
    X=x;
    Pxx=P;
}

void EKF::TimeUpdate(double time_,const Matrix & X_,const Matrix &STM)              //STM: State Tranform Matrix
{
    ++epoch;
    time=time_;
    X=X_;
    Pxx=STM*Pxx*Trans(STM);
}


void EKF::MeasUpdate(const Matrix &Z,const Matrix &g,const Matrix Weight_inv,const Matrix &G)        //Z:observations  //G:dz/dx
{   
    Matrix K;                                                                                              //g:model value by Propagating
    K=Pxx*Trans(G)*Inv_LU(Weight_inv+G*Pxx*Trans(G));
    X=X+K*(Z-g);
    Pxx=(Eye(n)-K*G)*Pxx;
}
