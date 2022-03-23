#pragma once 

#ifndef _GNSS_H
#define _GNSS_H
#include<cmath>
#include"VecMat.h"
struct Satllite
{
    unsigned int SatID;
    double XYZ_Sat[3];
    double parameter[6];
    /* data */
};

struct Site
{
    int SiteID;
    double XYZ_Site[3];
    /* data */
};


typedef struct BaseLine 
{
    int SiteID[2];
    double length;
    BaseLine()
    {
        length=0.0;
    }
    BaseLine(Site site1,Site site2)
    {
        SiteID[0]=site1.SiteID;
        SiteID[1]=site2.SiteID;
        length=sqrt(pow(site1.XYZ_Site[0]-site2.XYZ_Site[0],2)+pow(site1.XYZ_Site[1]-site2.XYZ_Site[1],2)+
               pow(site1.XYZ_Site[2]-site2.XYZ_Site[2],2));
    }
}BL;

struct DoubleAmb
{
    int site_id[2];
    int sat_id[2];
    int DDAmb;
    /* data */
};

class EPHbrdc
{
    public:
    EPHbrdc();
    ~EPHbrdc();
    void Init(const double* array,int n);
    void Init(std::string str);
    void CalPos(double gpstime,double* XY);

    void Pos2WGS84(double gpstime,double *XYZ);

    private:
    int _prn;
    //轨道根数
    double _gt;
    double _Omega;
    double _M0;
    double _e;
    double _SqrtA;
    double _omega;
    double _i;
    //轨道摄动参数
    double _DeltaN;
    double _DeltaOmega;
    double _DeltaI;
    double* _WaveU;
    double* _WaveI;
    double* _WaveR;
};


bool CheckIndepentBL(Matrix &IndependentMat,const BL &newline);

bool CheckIndepentAmb(Matrix &IndependentMat,const DoubleAmb & newAmb,int satnum,int sitenum);

void FindIndepBL(BL *inBL,BL* outBL,int total,int &independent,int Sitenum);

void FindIndepAmb(DoubleAmb *indoubleamb,DoubleAmb* outdoubleamb,int total,int &independentamb,int Sitenum,int Satnum);

int CompareBL(void *BL1,void *BL2);

void BLH2XYZ(const double* BLH,double* XYZ);


void XYZ2BLH(const double* XYZ,double* BLH);                      //CGCS-2000

Matrix TranMatofENU(const double* XYZ);                //********CGCS-2000**********//  

void EcceAno(double M,double e,double &b);
void TrueAno(double e,double E,double &f);

void BrdcReader(std::string filename, int &PRN,double* data);
void BrdcReader(std::string filename,int *PRN ,double** data);

void GetBrdcPara(double *data1,double *ephpara);

#endif