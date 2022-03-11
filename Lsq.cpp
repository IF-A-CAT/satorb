/*****************Least Square Adjustment with Weight Matrix**************
*this demo is for lsq with(without) weight matrix
*there is also a gross detection in the file,but only can detect one gross error
**********************************************/
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

using namespace std;

//*****************Least Square Adjustment **************//
int main()
{   
    fstream file_loader("../datafile/CUSV20210222_BDS.txt");
    ofstream file_out("../outputfile/GPSResult_XYZ_Weight_detection.txt");
    ofstream file_out1("../outputfile/GPSResult_Mean_Weight_detection.txt");
    file_out1<<"       X(m)     "<<"       Y(m)     "<<"         Y(m)         "<<endl; 
    std::string temp;
    int count=0;
    string info;
    double XYZ[3];
    double S0=0.0,Clk=0.0,obs[4],unit_error;
    Matrix x_before(4,1),x_post(4,1),Pxx(4,4),
            XYZ_Sat(3,1),XYZ_Site(3,1);
    x_post(1,0)=100; 
    // int count_n=0;                                       //init for compare
    while(getline(file_loader,temp))
    {
        count++;
        int num;
        std::string infomation=temp;
        file_out<<temp<<endl;
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
        Matrix WeightMat(num,num);
        Matrix l(num,1),B(num,4),V,L(num,num);
        for(int i=0;i<num;i++)                                        //get B and l
        {   
            getline(file_loader,temp);

            /************this part is for weight matrix*********/
            //cout<<static_cast<int>(temp[2])<<endl;
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
            /*****************************************************/
            Split(temp,' ',obs);
            XYZ_Sat(0,0)=obs[0];
            XYZ_Sat(1,0)=obs[1];
            XYZ_Sat(2,0)=obs[2];
            
            ObsMat(i,0)=obs[0];ObsMat(i,1)=obs[1];ObsMat(i,2)=obs[2];ObsMat(i,3)=obs[3];
            S0=sqrt(pow(XYZ_Site(0,0)-XYZ_Sat(0,0),2)+pow(XYZ_Site(1,0)-XYZ_Sat(1,0),2)+pow(XYZ_Site(2,0)-XYZ_Sat(2,0),2))+Clk;
            l(i,0)=ObsMat(i,3)-S0;
            B(i,0)=(XYZ_Site(0,0)-XYZ_Sat(0,0))/S0;B(i,1)=(XYZ_Site(1,0)-XYZ_Sat(1,0))/S0;B(i,2)=(XYZ_Site(2,0)-XYZ_Sat(2,0))/S0;B(i,3)=1.0;
    
        }
        //add weight
        B=L*B;
        l=L*l;
        int n=1;
        double delta_x=100;
        while(1)
        {
            n++;
            if(n>20||delta_x<0.1)
            {
                break;
            }
            // cout<<l<<B<<endl;
            Lsq lsq_fliter(l,B,4);
            XYZ_Site=XYZ_Site+lsq_fliter.get_x().slice(0,2,0,0);
            // cout<<lsq_fliter.x.slice(0,2,0,0);
            unit_error=lsq_fliter.get_corr();
            Aux=lsq_fliter.get_pxx();
            V=lsq_fliter.get_v();
            delta_x=sqrt(Norm(lsq_fliter.get_x()));
            Clk+=lsq_fliter.get_x()(3,0);
            for(int i=0;i<num;i++)
            {
                S0=sqrt(pow(XYZ_Site(0,0)-ObsMat(i,0),2)+pow(XYZ_Site(1,0)-ObsMat(i,1),2)+pow(XYZ_Site(2,0)-ObsMat(i,2),2))+Clk;
                l(i,0)=ObsMat(i,3)-S0;
                B(i,0)=(XYZ_Site(0,0)-ObsMat(i,0))/S0;B(i,1)=(XYZ_Site(1,0)-ObsMat(i,1))/S0;B(i,2)=(XYZ_Site(2,0)-ObsMat(i,2))/S0;
            }
            // add weight
            B=L*B;
            l=L*l;
        }




        //***********this part is for gross_error_detection and remove gross*************
        if(gross_error_detection(num-4,1.0,V,Eye(num)))
        {   
            for(int i=0;i<num;i++)
            {   
                Matrix Aux2(ObsMat);
                Aux2.exchange_row(i,num-1);
                Matrix Aux3(L);
                Aux3.exchange_row(i,num-1);
                Matrix H(num-1,4);
                Matrix Z(num-1,1);
                for(int i=0;i<5;i++)
                {
                Matrix Aux1=Aux2.slice(0,num-2,0,3);
                get_coef_mat(Aux1,XYZ_Site,Clk,H,Z);
                // cout<<H<<Z<<endl;
                Lsq least_square(Aux3.slice(0,num-2,0,num-2)*Z,Aux3.slice(0,num-2,0,num-2)*H,4);
                Clk+=least_square.get_x()(3,0);
                Aux=least_square.get_pxx();
                unit_error=least_square.get_corr();
                V=least_square.get_v();
                XYZ_Site=XYZ_Site+least_square.get_x().slice(0,2,0,0);
                }
                // cout<<V<<endl;
                if(!gross_error_detection(num-5,1.0,V,Eye(num-1)))
                {
                    
                    cout<<infomation<<endl;
                    cout<<"Gross_error is at the "<<i+1<<"th Observation"<<endl<<endl;
                    break;
                }
            }
        }
        /*****************************************************/

        //file output
        file_out<<"       X(m)       "<<"       Y(m)       "<<"     Z(m)     "<<"         T(m)      "<<endl;
        file_out<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
        fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
        fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<
        fixed<<setprecision(4)<<setw(16)<<Clk<<endl;
        file_out1<<fixed<<setprecision(4)<<setw(16)<<XYZ_Site(0,0)<<
        fixed<<setprecision(4)<<setw(16)<<XYZ_Site(1,0)<<
        fixed<<setprecision(4)<<setw(16)<<XYZ_Site(2,0)<<endl;
        file_out<<"验后单位权中误差:  "<<fixed<<setprecision(4)<<setw(16)<<unit_error<<endl;
        file_out<<"验后估计方差:"<<endl;
        file_out<<Aux*pow(unit_error,2);
        file_out<<"------------------------------------------------------------------------------------------------"<<endl;
}
}