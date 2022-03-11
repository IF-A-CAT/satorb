#include<iostream>
#include<fstream>
#include<cmath>
#include"./include/VecMat.h"
#include"./include/Algorithm.h"
#include<string>
#include"./include/Const.h"
#include"./include/GNSS.h"
#include<iomanip>
#include<ctime>
// #include<io.h>
// #include<direct.h>
using namespace std;




//****************Square Root Information Filter******************//
int main()
{   clock_t start,end;
    start=clock();
    fstream outfile("../GPSResult_Srif.txt",std::ios::out);                 //file will be created only without instream mode
    std::fstream fileread("../CUSV20210222_BDS.txt");
    std::string temp;
    int count=0;
    string info;
    double S0=0.0,Clk=0.0,obs[4],unit_error;
    Matrix x_before(4,1),x_post(4,1),Pxx(4,4),
            XYZ_Sat(3,1),XYZ_Site(3,1);
    x_post(1,0)=100;                                        //init for compare
    while(getline(fileread,temp))
    {
        count++;
        int num;
        outfile<<temp<<endl;
        for(int len=temp.length()-1;len>=0;len--)
        {
            if(temp[len]!=' '&&temp[len]!='\0'&&temp[len]!='\n')
            {
                // cout<<temp.substr(len-2,3)<<endl;
                num=atoi(temp.substr(len-3,3).c_str());
                break;
            }
        }                               //print information
        Matrix Aux,ObsMat(num,4);
        Matrix L(num,num);
        Matrix l(num,1),B(num,4);
        Clk=0.0;
        for(int i=0;i<num;i++)                                        //get B and l
        {   
            getline(fileread,temp);
            if((temp[1]=='0'&&int(temp[2])-48<=5)||atoi(temp.substr(1,2).c_str())==59||temp[1]=='6')
                {
                    L(i,i)=1.0/3.0;
                }
                else if(atoi(temp.substr(1,2).c_str())==31||atoi(temp.substr(1,2).c_str())==56||atoi(temp.substr(1,2).c_str())==57||
                atoi(temp.substr(1,2).c_str())==58||atoi(temp.substr(1,2).c_str())==38||atoi(temp.substr(1,2).c_str())==39||atoi(temp.substr(1,2).c_str())==40||
                atoi(temp.substr(1,2).c_str())==13||atoi(temp.substr(1,2).c_str())==16||(temp[1]=='0'&&int(temp[2])-48>5)||atoi(temp.substr(1,2).c_str())==10)
                {
                    L(i,i)=1.0/2.0;
                }
                else{
                    L(i,i)=1.0;
                }
            Split(temp,' ',obs);
            XYZ_Sat(0,0)=obs[0];
            XYZ_Sat(1,0)=obs[1];
            XYZ_Sat(2,0)=obs[2];
            // L(i,0)=obs[3];
            ObsMat(i,0)=obs[0];ObsMat(i,1)=obs[1];ObsMat(i,2)=obs[2];ObsMat(i,3)=obs[3];
            S0=sqrt(pow(XYZ_Site(0,0)-XYZ_Sat(0,0),2)+pow(XYZ_Site(1,0)-XYZ_Sat(1,0),2)+pow(XYZ_Site(2,0)-XYZ_Sat(2,0),2))+Clk;
            l(i,0)=ObsMat(i,3)-S0;
            B(i,0)=(XYZ_Site(0,0)-XYZ_Sat(0,0))/S0;B(i,1)=(XYZ_Site(1,0)-XYZ_Sat(1,0))/S0;B(i,2)=(XYZ_Site(2,0)-XYZ_Sat(2,0))/S0;B(i,3)=1.0;
        }
            int n=0;
            B=L*B;
            l=L*l;
            while(1)                            //Recursion
            {
                n++;
                if(count==1)
                {   
                    Lsq lsq(l,B,4);
                    // x_before=x_post;

                    x_post=lsq.get_x();
                    unit_error=lsq.get_corr();
                    Aux=lsq.get_pxx();
                    // cout<<x_post-Inv_LU(Trans(B)*B)*Trans(B)*l;
                    XYZ_Site=XYZ_Site+lsq.get_x().slice(0,2,0,0);
                    Clk+=lsq.get_x()(3,0);
                    // Update B and l
                    for(int i=0;i<num;i++)
                    {
                        S0=sqrt(pow(XYZ_Site(0,0)-ObsMat(i,0),2)+pow(XYZ_Site(1,0)-ObsMat(i,1),2)+pow(XYZ_Site(2,0)-ObsMat(i,2),2))+Clk;
                        l(i,0)=ObsMat(i,3)-S0;
                        B(i,0)=(XYZ_Site(0,0)-ObsMat(i,0))/S0;B(i,1)=(XYZ_Site(1,0)-ObsMat(i,1))/S0;B(i,2)=(XYZ_Site(2,0)-ObsMat(i,2))/S0;
                    }
                    
                }
                else
                {
                    Lsq lsq(Pxx,l,B,4);
                    // x_before=x_post;
                    x_post=lsq.get_x();
                    unit_error=lsq.get_corr();
                    Aux=lsq.get_pxx();
                    XYZ_Site=XYZ_Site+lsq.get_x().slice(0,2,0,0);
                    Clk+=lsq.get_x()(3,0);
                     // Update B and l
                    for(int i=0;i<num;i++)
                    {
                        S0=sqrt(pow(XYZ_Site(0,0)-ObsMat(i,0),2)+pow(XYZ_Site(1,0)-ObsMat(i,1),2)+pow(XYZ_Site(2,0)-ObsMat(i,2),2))+Clk;
                        l(i,0)=ObsMat(i,3)-S0;
                        B(i,0)=(XYZ_Site(0,0)-ObsMat(i,0))/S0;B(i,1)=(XYZ_Site(1,0)-ObsMat(i,1))/S0;B(i,2)=(XYZ_Site(2,0)-ObsMat(i,2))/S0;
                    }
                    // Pxx=lsq.Pxx_post;
                }
                // cout<<l<<'\n'<<B<<endl;
                // cout<<<<endl;Inv_LU(Pxx_)
                if(Norm(x_post)<0.1||n>15)          //迭代20次限制与完全达到1e-8要求之间差异均小于1mm，且加入20次限制可以明显提高cpu计算速度，减少时间
                {
                    Pxx=Aux;
                    break;
                }
            }                                                  //problem data: 16   1086    
            outfile<<"       X(m)       "<<"       Y(m)       "<<"     Z(m)     "<<"         T(m)      "<<endl;
            outfile<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
                  fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
                  fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<
                  fixed<<setprecision(4)<<setw(16)<<Clk<<endl;
            outfile<<"验后单位权中误差:  "<<fixed<<setprecision(4)<<setw(16)<<unit_error<<endl;
            outfile<<"验后估计方差:"<<endl;
            outfile<<Pxx*unit_error*unit_error;
            outfile<<"------------------------------------------------------------------------------------------------"<<endl;
            // cout<<XYZ_Site<<'\n'<<Clk<<endl;
            Pxx(3,3)=1.0;
            for(int i=0;i<3;i++)
            {
                Pxx(i,3)=0.0;
            }
            for(int j=0;j<3;j++)
            {
                Pxx(3,j)=0.0;
            }
        }
    end=clock();
    outfile<<"Total run time :  "<<double(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
//****************end***************//