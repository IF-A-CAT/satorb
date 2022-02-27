#include<iostream>
#include<fstream>
#include<cmath>
#include"VecMat.h"
#include"Algorithm.h"
#include<string>
#include"Const.h"
#include"GNSS.h"
#include<iomanip>
#include<ctime>
// #include<io.h>
// #include<direct.h>
using namespace std;



// int main()
// {   
//     fstream file_loader("E:\\课程学习文件\\最优估计\\data_adjustment\\CUSV20210222_BDS.txt");
//     ofstream file_out("E:\\课程学习文件\\最优估计\\data_adjustment\\GPSResult_test.txt");
//     ofstream file_out1("E:\\课程学习文件\\最优估计\\data_adjustment\\GPSResult_test1.txt");
//     file_out1<<"       X(m)     "<<"       Y(m)     "<<"         Y(m)         "<<endl; 
//     std::string temp;
//     int count=0;
//     string info;
//     double XYZ[3];
//     double S0=0.0,Clk=0.0,obs[4],unit_error;
//     Matrix x_before(4,1),x_post(4,1),Pxx(4,4),
//             XYZ_Sat(3,1),XYZ_Site(3,1);
//     x_post(1,0)=100;                                        //init for compare
//     while(getline(file_loader,temp))
//     {
//         count++;
//         int num;
//         file_out<<temp<<endl;
//         for(int len=temp.length()-1;len>=0;len--)
//         {
//             if(temp[len]!=' '&&temp[len]!='\0'&&temp[len]!='\n')
//             {
//                 // cout<<temp.substr(len-2,3)<<endl;
//                 num=atoi(temp.substr(len-2,3).c_str());
//                 break;
//             }
//         }                               //print information
//         Matrix Aux,ObsMat(num,4);
//         Matrix L(num,num);
//         Matrix l(num,1),B(num,4);
//         for(int i=0;i<num;i++)                                        //get B and l
//         {   
//             getline(file_loader,temp);
//             // cout<<static_cast<int>(temp[2])<<endl;
//             if((temp[1]=='0'&&int(temp[2])-48<=5)||atoi(temp.substr(1,2).c_str())==59||temp[1]=='6')
//             {
//                 L(i,i)=1.0/3.0;
//             }
//             else if(atoi(temp.substr(1,2).c_str())==31||atoi(temp.substr(1,2).c_str())==56||atoi(temp.substr(1,2).c_str())==57||
//             atoi(temp.substr(1,2).c_str())==58||atoi(temp.substr(1,2).c_str())==38||atoi(temp.substr(1,2).c_str())==39||atoi(temp.substr(1,2).c_str())==40||
//             atoi(temp.substr(1,2).c_str())==13||atoi(temp.substr(1,2).c_str())==16||(temp[1]=='0'&&int(temp[2])-48>5)||atoi(temp.substr(1,2).c_str())==10)
//             {
//                 L(i,i)=1.0/2.0;
//             }
//             else{
//                 L(i,i)=1.0;
//             }
//             Split(temp,' ',obs);
//             XYZ_Sat(0,0)=obs[0];
//             XYZ_Sat(1,0)=obs[1];
//             XYZ_Sat(2,0)=obs[2];
//             // L(i,0)=obs[3];
//             ObsMat(i,0)=obs[0];ObsMat(i,1)=obs[1];ObsMat(i,2)=obs[2];ObsMat(i,3)=obs[3];
//             S0=sqrt(pow(XYZ_Site(0,0)-XYZ_Sat(0,0),2)+pow(XYZ_Site(1,0)-XYZ_Sat(1,0),2)+pow(XYZ_Site(2,0)-XYZ_Sat(2,0),2))+Clk;
//             l(i,0)=ObsMat(i,3)-S0;
//             B(i,0)=(XYZ_Site(0,0)-XYZ_Sat(0,0))/S0;B(i,1)=(XYZ_Site(1,0)-XYZ_Sat(1,0))/S0;B(i,2)=(XYZ_Site(2,0)-XYZ_Sat(2,0))/S0;B(i,3)=1.0;
//         }
//         int n=1;
//         B=L*B;
//         l=L*l;
//         double delta_x=100;
//         while(1)
//         {
//             if(n>15||delta_x<0.1)
//             {
//                 break;
//             }
//             n++;
//             Lsq lsq_fliter(l,B,4);
//             XYZ_Site=XYZ_Site+lsq_fliter.get_x().slice(0,2,0,0);
//             // cout<<lsq_fliter.x.slice(0,2,0,0);
//             unit_error=lsq_fliter.get_corr();
//             Aux=lsq_fliter.get_pxx();
//             delta_x=sqrt(Norm(lsq_fliter.get_x()));
//             Clk+=lsq_fliter.get_x()(3,0);
//             for(int i=0;i<num;i++)
//             {
//                 S0=sqrt(pow(XYZ_Site(0,0)-ObsMat(i,0),2)+pow(XYZ_Site(1,0)-ObsMat(i,1),2)+pow(XYZ_Site(2,0)-ObsMat(i,2),2))+Clk;
//                 l(i,0)=ObsMat(i,3)-S0;
//                 B(i,0)=(XYZ_Site(0,0)-ObsMat(i,0))/S0;B(i,1)=(XYZ_Site(1,0)-ObsMat(i,1))/S0;B(i,2)=(XYZ_Site(2,0)-ObsMat(i,2))/S0;
//             }

//         }
            
//         // file_out<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
//         // fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
//         // fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<endl;
//         // XYZ[0]+=XYZ_Site(0,0);
//         // XYZ[1]+=XYZ_Site(1,0);
//         // XYZ[2]+=XYZ_Site(2,0);
//         file_out<<"         X(m)      "<<"           Y(m)         "<<"         Z(m)        "<<"          T(m)       "<<endl;
//         file_out<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
//         fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
//         fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<
//         fixed<<setprecision(8)<<setw(16)<<Clk<<endl;
//         file_out1<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
//         fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
//         fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<endl;
//         file_out<<"验后单位权中误差:  "<<fixed<<setprecision(4)<<setw(16)<<unit_error<<endl;
//         file_out<<"验后估计方差:"<<endl;
//         file_out<<Aux*pow(unit_error,2);
//         file_out<<"------------------------------------------------------------------------------------------------"<<endl;
// }
        // file_out<<"-----------------Mean XYZ(m)------------------"<<endl;
        // file_out<<fixed<<setprecision(4)<<setw(16)<<XYZ[0]/double(count)<<
        // fixed<<setprecision(4)<<setw(16)<<XYZ[1]/double(count)<<
        // fixed<<setprecision(4)<<setw(16)<<XYZ[2]/double(count)<<endl;

// }
    // double array[3]={4,2,1};
    // int i=0;
    // cout<<array[i++]<<" "<<array[++i];
// {
//     BL baseline[10];
//     for(auto &i:baseline)
//     {
//         i.length=rand();
//         cout<<i.length<<" ";
//     }
//     BL outbaseline[10];
//     int total=10,temp=0;
//     FindIndepBL(baseline,outbaseline,total,temp);                                                             //Check qsort
//     for(auto &i:baseline)
//     {
//         // i.length=rand();
//         cout<<i.length<<" ";
//     }
// }
// {
//     double array[100]={92,  99   ,  1 ,    8,    15,    67,    74,    51,    58,    40,
//                       98,  80   ,  7 ,   14,    16,    73,    55,    57,    64,    41,
//                        4,  81   , 88 ,   20,    22,    54,    56,    63,    70,    47,
//                       85,  87   , 19 ,   21,     3,    60,    62,    69,    71,    28,
//                       86,  93   , 25 ,    2,     9,    61,    68,    75,    52,    34,                       //check QR
//                       17,  24   , 76 ,   83,    90,    42,    49,    26,    33,    65,
//                       23,   5   , 82 ,   89,    91,    48,    30,    32,    39,    66,
//                       79,   6   , 13 ,   95,    97,    29,    31,    38,    45,    72,
//                       10,  12   , 94 ,   96,    78,    35,    37,    44,    46,    53,
//                       11,  18   ,100 ,   77,    84,    36,    43,    50,    27,    59};
//     double array1[4]={1,1,1,1};
//     Matrix l(array1,4,1);
//     // Matrix M1(M|l);
//     Matrix M(array,10,10);
//     if(IsColFullRank(M))
//     {
//         cout<<"yes"<<endl;
//     }
//     cout<<Inv_LU(M)<<endl;
// }

// int main()
// {
//     double array[16]={0.9177 ,   -1.0404,    -0.9598 ,   -0.9903,
//                     -1.0404 ,    2.7034   ,  1.2022   ,  1.8405,
//                     -0.9598  ,   1.2022   ,  1.9910  ,   1.3921,
//                     -0.9903    , 1.8405    , 1.3921 ,    1.6431};
//    Matrix A(array,4,4);
//    Matrix *M;
//    M=QR_decompose(A);
//    cout<<M[0]<<'\n'<<M[1]<<endl;
//    cout<<M[0]*M[1]<<endl;
//    cout<<Inv(M[1])<<endl;

//    return 0;
// }
//****************Root mean square information filtering***************//
int main()
{   clock_t start,end;
    start=clock();
    fstream outfile("E:\\课程学习文件\\最优估计\\data_adjustment\\GPSResult_Info.txt",std::ios::out);                 //file will be created only without instream mode
    std::fstream fileread("E:\\课程学习文件\\最优估计\\data_adjustment\\CUSV20210222_BDS.txt");
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
                num=atoi(temp.substr(len-2,3).c_str());
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
                // if(n>10)                            //迭代20次限制与完全达到1e-8要求之间差异均小于1mm，且加入20次限制可以明显提高cpu计算速度，减少时间
                // {
                //     cerr<<"ERROR:Convergence failed!!!"<<endl;
                //     break;
                // }
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
                if(Norm(x_post)<0.1||n>15)
                {
                    Pxx=Aux;
                    break;
                }
            }                                                   //problem data: 16   1086    
            outfile<<"          X(m)         "<<"         Y(m)         "<<"         Z(m)       "<<"         T(m)         "<<endl;
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