#include"../include/GNSS.h"
#include"../include/VecMat.h"
#include<iostream>
#include<vector>
#include<string>

bool CheckIndepentBL(Matrix &IndependentMat,const BL &newline)                    //Site Index from 1
{
    Matrix column(IndependentMat.get_row(),1);
    column(newline.SiteID[0]-1,0)=1;
    column(newline.SiteID[1]-1,0)=-1;
    Matrix M(IndependentMat|column);
    if(IsColFullRank(M))
    {
        IndependentMat=M;
        return true;
    }
    else
        return false;
}


bool CheckIndepentAmb(Matrix &IndependentMat,const DoubleAmb & newAmb,int satnum,int sitenum)
{
    Matrix column(satnum*sitenum,1);
    column((newAmb.site_id[0]-1)*satnum+newAmb.sat_id[0]-1,0)=1;
    column((newAmb.site_id[1]-1)*satnum+newAmb.sat_id[0]-1,0)=-1;
    column((newAmb.site_id[0]-1)*satnum+newAmb.sat_id[1]-1,0)=-1;
    column((newAmb.site_id[1]-1)*satnum+newAmb.sat_id[1]-1,0)=1;
    Matrix M(IndependentMat|column);
    if(IsColFullRank(M))
    {
        IndependentMat=M;
        return true;
    }
    else
        return false;
}

int CompareBL(const void *BL1,const void *BL2)                      //To sort BL by length
{
    return int(((BL*)BL1)->length-((BL*)BL2)->length);
}

void FindIndepBL(BL *inBL,BL* outBL,int total,int &independentnum,int Sitenum)        //Find independentline use shared_ptr or delete outBL by yourself
{
    qsort(inBL,total,sizeof(BL),&CompareBL);                                 //sort inBL
    Matrix checkmat(Sitenum,2);
    // Matrix Aux;
    independentnum=0;
    std::vector<BL> outbl;                                            //Matrix used to check independence by column rank
    for(int i=0;i<total;i++)
    {
        if(i<2)
        {
            checkmat(inBL[i].SiteID[0]-1,independentnum++)=1;
            checkmat(inBL[i].SiteID[1]-1,independentnum++)=-1;
            outbl.push_back(inBL[i]);        //input two shortest BL
        }
        else
        {
            //--------------another way--------------------
            // Matrix AddCol(Sitenum,1);                                       //new col added to checkmatrix
            // AddCol(inBL[i].SiteID[0]-1,0)=1;
            // AddCol(inBL[i].SiteID[1]-1,0)=-1;
            // Aux= checkmat|AddCol;
            // if(IsColFullRank(Aux))                                          //check if the new BL is independent
            // {
            //     checkmat=Aux;
            //     ++independentnum;
            //     outbl.push_back(inBL[i]);
            // }
            //--------------------------------------------------
            if(CheckIndepentBL(checkmat,inBL[i]))
            {
                ++independentnum;
                outbl.push_back(inBL[i]);
            }
        }
    }
    delete [] outBL;
    outBL=new BL [independentnum];
    for(int i=0;i<independentnum;++i)
    {
        outBL[i].length=outbl[i].length;
        outBL[i].SiteID[0]=outbl[i].SiteID[0];
        outBL[i].SiteID[1]=outbl[i].SiteID[1];
    }
}


void FindIndepAmb(DoubleAmb *inDoubleAmb,DoubleAmb* outDoubleAmb,int total,int &independentnum,int Sitenum,int Satnum)
{   
    // qsort(indoubleamb,total,sizeof(DoubleAmb),&CompareAmb);                           //sort by the weight of each Amb
    Matrix Aux;
    independentnum=0;
    std::vector<DoubleAmb> outDA; 
    Matrix checkmat(Sitenum*Satnum,2);                                           //Matrix used to check independence by column rank
    for(int i=0;i<total;i++)
    {
        if(i<2)
        {
            // Matrix column(Satnum*Sitenum,1);
            checkmat((inDoubleAmb[i].site_id[0]-1)*Satnum+inDoubleAmb[i].sat_id[0]-1,independentnum++)=1;
            checkmat((inDoubleAmb[i].site_id[1]-1)*Satnum+inDoubleAmb[i].sat_id[0]-1,independentnum++)=-1;
            checkmat((inDoubleAmb[i].site_id[0]-1)*Satnum+inDoubleAmb[i].sat_id[1]-1,independentnum++)=-1;
            checkmat((inDoubleAmb[i].site_id[1]-1)*Satnum+inDoubleAmb[i].sat_id[1]-1,independentnum++)=1;   
            outDA.push_back(inDoubleAmb[i]);        //input two shortest BL
        }
        else
        {
            if(CheckIndepentAmb(checkmat,inDoubleAmb[i],Satnum,Sitenum))                                          //check if the new BL is independent
            {
                ++independentnum;
                outDA.push_back(inDoubleAmb[i]);
            }
        }
    }
    delete [] outDoubleAmb;                                         //when recursion,prevent 
    outDoubleAmb=new DoubleAmb [independentnum];
    for(int i=0;i<independentnum;++i)
    {
        outDoubleAmb[i].DDAmb=outDA[i].DDAmb;
        outDoubleAmb[i].sat_id[0]=outDA[i].sat_id[0];
        outDoubleAmb[i].sat_id[1]=outDA[i].sat_id[1];
        outDoubleAmb[i].site_id[0]=outDA[i].site_id[0];
        outDoubleAmb[i].site_id[1]=outDA[i].site_id[1];
    }
}
