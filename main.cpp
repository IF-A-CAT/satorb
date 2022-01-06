#include<iostream>
#include"VecMat.h"
using namespace std;

int main(void)
{
    Vector V1(2,1,3);
    Vector V2(3,2,1);
    Vector V3(3);
    double array[20]={1,2,3,4,5,6,7,8,9,10,1,3,4,2,3,6,4,5,3,10};
    Matrix M1(array,4,5);
    Matrix M2(4,4);
    Matrix* M3;
    M3=QR_decompose(M1);

    cout<<M3[0]<<M3[1]<<M3[0]*Trans(M3[0]);
    cout<<M3[0]*M3[1];
    return 0;
}