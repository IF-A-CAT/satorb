#include<iostream>
#include<cmath>
#include"VecMat.h"
using namespace std;

int main()
{   
    Vector V1(2,1,3);
    Vector V2(3,2,1);
    Vector V3(3);
    double array[16]={1,2,3,4,5,6,7,8,9,10,1,3,4,2,3,5};
    Matrix M1(array,4,4);
    double pi=M_PI;
    Matrix Rx;
    Rx=Full(2.0,1);
    cout<<Rx;
    // Matrix M2(4,4);
    // Matrix M3;
    // Rx(M_PI);
    // M3=Inv_LU(M1);
    // M2=Inv_QR(M1);
    // cout<<M3<<"\n"<<M2<<M2-M3;
    // // cout<<M3[0]<<M3[1]<<M3[0]*Trans(M3[0]);
    // // cout<<M3[0]*M3[1];
    return 0;
}