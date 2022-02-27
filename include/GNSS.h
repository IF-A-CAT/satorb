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


bool CheckIndepentBL(Matrix &IndependentMat,const BL &newline);

bool CheckIndepentAmb(Matrix &IndependentMat,const DoubleAmb & newAmb,int satnum,int sitenum);

void FindIndepBL(BL *inBL,BL* outBL,int total,int &independent,int Sitenum);

void FindIndepAmb(DoubleAmb *indoubleamb,DoubleAmb* outdoubleamb,int total,int &independentamb,int Sitenum,int Satnum);

int CompareBL(void *BL1,void *BL2);
#endif