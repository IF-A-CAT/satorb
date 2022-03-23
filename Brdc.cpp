#include"./include/GNSS.h"
// #include"./include/VecMat.h"
// #include"./include/Const.h"
#include<cmath>
#include<iomanip>
#include<fstream>

int main()
{
    std::string filename="./datafile/test_eph";
    std::ofstream outputfile("./outputfile/EPH.txt");
    EPHbrdc eph;
    double XYZ[3]={0.0};
    eph.Init(filename);
    eph.Pos2WGS84(79439.9315768461,XYZ);
    outputfile<<" Epoch "<<"    X(m)      "<<"      Y(m)      "<<"      Z(m)      "<<"\n";
    outputfile<<"  "<<1<<"  ";
    for(int i=0;i<3;i++)
    {
        outputfile<<std::fixed<<std::setw(16)<<std::setprecision(6)<<XYZ[i];
    }
    outputfile<<"\n";
    outputfile<<"  "<<2<<"  ";
    eph.Pos2WGS84(86279.9237728612,XYZ);
    for(int i=0;i<3;i++)
    {
        outputfile<<std::fixed<<std::setw(16)<<std::setprecision(6)<<XYZ[i];
    }
    outputfile<<"\n";
    outputfile<<"  "<<3<<"  ";
    eph.Pos2WGS84(83099.9292294162,XYZ);
    for(int i=0;i<3;i++)
    {
        outputfile<<std::fixed<<std::setw(16)<<std::setprecision(6)<<XYZ[i];
    }
    return 0;
}