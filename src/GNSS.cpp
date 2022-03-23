#include"../include/GNSS.h"
#include"../include/Const.h"
#include"../include/VecMat.h"
#include<iostream>
#include<fstream>
#include<vector>
#include<string.h>

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

void BLH2XYZ(const double* BLH,double* XYZ)
{
    double a=6378137;
    double b=6356752.31414;
    double N,f,e;
    e=sqrt(a*a-b*b)/a;
    N=a/sqrt(1-e*e*sin(BLH[0])*sin(BLH[0]));
    XYZ[0]=(N+BLH[2])*cos(BLH[0])*cos(BLH[1]);
    XYZ[1]=(N+BLH[2])*cos(BLH[0])*sin(BLH[1]);
    XYZ[2]=(N*(1-e*e)+BLH[2])*sin(BLH[0]);
}

void XYZ2BLH(const double* XYZ,double* BLH)                       
{
    double a=6378137;
    double b= 6356752.31414;
    double e=sqrt(a*a-b*b)/a;
    double TanB=0,Btemp=1;
    while(fabs(Btemp-TanB)>1e-6)
    {
        Btemp=TanB;
        TanB=(1/sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1]))*(XYZ[2]+(a*e*e*TanB)/sqrt(1+TanB*TanB-e*e*TanB*TanB));
    }
    BLH[1]=atan(XYZ[1]/XYZ[0]);
    if(XYZ[0]>=0)
    {
        if(XYZ[1]<0)
        {
            BLH[1]=2*pi+BLH[1];
        }
    }
    else if(XYZ[0]<0)
    {
            BLH[1]=BLH[1]+pi;
    }
    BLH[0]=atan(TanB);
    BLH[2]=sqrt(XYZ[0]*XYZ[0]+XYZ[1]*XYZ[1])/cos(BLH[0])-a/sqrt(1-e*e*sin(BLH[0])*sin(BLH[0]));
}
Matrix TranMatofENU(const double* XYZ)                   //********CGCS-2000**********//  
{
    Matrix Tran(3,3);
    double BLH[3];
    XYZ2BLH(XYZ,BLH);
    Tran=R_z(BLH[1]+pi/2);
    Tran=R_x(pi/2-BLH[0])*Tran;
    return Tran;
}

EPHbrdc::EPHbrdc()
{
    _WaveI=new double [2];
    _WaveR=new double [2];
    _WaveU=new double [2];
}
EPHbrdc::~EPHbrdc()
{
    delete [] _WaveI;
    delete [] _WaveU;
    delete [] _WaveR;
}
void EPHbrdc::Init(const double* array,int n)
{
    if(n==6)
    {
        _M0=array[0];
        _e=array[1];
        _SqrtA=array[2];
        _i=array[3];
        _omega=array[4];
        _Omega=array[5];
    }
    else if (n==16)
    {
        _M0=array[0];
        _e=array[1];
        _SqrtA=array[2];
        _i=array[3];
        _omega=array[4];
        _Omega=array[5];
        _DeltaN=array[6];
        _DeltaOmega=array[7];
        _DeltaI=array[8];
        _WaveU[0]=array[9];
        _WaveU[1]=array[10];
        _WaveI[0]=array[11];
        _WaveI[1]=array[12];
        _WaveR[0]=array[13];
        _WaveR[1]=array[14];
        _gt=array[15];
    }
    else{
        std::cerr<<"ERROR:EPHbrdc INPUT fail ! \n";
        exit(1);
    }
}
//newton iteration
void EcceAno(double M,double e,double &E)
{
    double aux,temp,deriv,deltaE=1.0;
    E=M;
    while(fabs(deltaE)>1e-12)
    {
    aux=E;
    temp=E-e*sin(E)-M;
    deriv=1-e*cos(E);
    E=E-temp/deriv;
    deltaE=aux-E;
    }
}

void TrueAno(double e,double E,double &f)
{
    double temp1,temp2;
    temp1=cos(E)-e;
    temp2=sqrt(1-e*e)*sin(E);
    if(temp1>0&&temp2>=0)
    {
        f=atan(temp2/temp1);
    }
    else if(temp1<0)
    {
        f=atan(temp2/temp1)+pi;
    }
    else if(temp1>0&&temp2<0)
    {
        f=atan(temp2/temp1)+2*pi;
    }
}

void EPHbrdc::CalPos(double gpstime,double* XY)
{
    double M,a,E,f,u,r;
    double delta_t=gpstime-_gt;
   //get f/u
    M=_M0+(sqrt(GM)/pow(_SqrtA,3)+_DeltaN)*(gpstime-_gt);
    EcceAno(M,_e,E);
    TrueAno(_e,E,f);
    
    u=_omega+f;
    
    
    //get r
    r=_SqrtA*_SqrtA*(1-_e*cos(E));
    
    //correction r
    r=r+_WaveR[0]*cos(2*u)+_WaveR[1]*sin(2*u);

    //calculate XY
    //correction u
    u=u+_WaveU[0]*cos(2*u)+_WaveU[1]*sin(2*u);
    XY[0]=r*cos(u);
    XY[1]=r*sin(u);
}

void EPHbrdc::Pos2WGS84(double gpstime,double *XYZ)
{
    double M,i,E,f,L,PosXY[2],u;
    CalPos(gpstime,PosXY);
    M=_M0+(sqrt(GM)/pow(_SqrtA,3)+_DeltaN)*(gpstime-_gt);
    EcceAno(M,_e,E);
    TrueAno(_e,E,f);
    u=_omega+f;
    
    //correction i
    i = _i + _DeltaI * ( gpstime - _gt ) + _WaveI[0] * cos(2*u) + _WaveI[1] * sin(2*u);
    //calculate longtitute
    L=_Omega+(_DeltaOmega-w) * gpstime -_DeltaOmega * _gt ;

    //correction u
    // u=u+_WaveU[0]*cos(2*u)+_WaveU[1]*sin(2*u);
    //transform matrix
    Matrix rotate,Pos(3,1),WGSPos;
    rotate=R_x(-i);
    rotate= R_z(-L) * rotate ;
    // std::cout<<rotate;
    //position of satellite in plane of satellite
    Pos(0,0)=PosXY[0];
    Pos(1,0)=PosXY[1];
    Pos(2,0)=0.0;
    
    //WGS-84 location 
    WGSPos=rotate * Pos;

    // std::cout<<WGSPos;
    XYZ[0]=WGSPos(0,0);
    XYZ[1]=WGSPos(1,0);
    XYZ[2]=WGSPos(2,0);
}

void EPHbrdc::Init(std::string str)
{
    int prn;
    double data[31],ephpara[16];
    BrdcReader(str,prn,data);
    GetBrdcPara(data,ephpara);
    Init(ephpara,16);
    _prn=prn;
}

void BrdcReader(std::string filename, int &PRN,double* data)
{
    std::string str;
    int leapsec;
    std::ifstream loader(filename,std::ios::in);
    //HEADER READER
    while(getline(loader,str))
    {
        const char* str_copy=str.c_str();
        if( NULL != strstr(str_copy,"LEAP SECONDS"))
        {
            leapsec = atoi(str.substr(0,60).c_str()); 
        }
        else if(NULL != strstr(str_copy,"END OF HEADER"))
        {
            break;
        }
    }
    //DATA READER
    while(getline(loader,str))
    {
        PRN=atoi(str.substr(0,3).c_str());
        int year,month,day,hour,min;
        float sec;
        // double tempdata[4]={0.0};

        year=atoi(str.substr(3,2).c_str());
        month=atoi(str.substr(6,2).c_str());
        day=atoi(str.substr(9,2).c_str());
        hour=atoi(str.substr(12,2).c_str());
        min=atoi(str.substr(15,2).c_str());
        sec=atof(str.substr(18,4).c_str());
/**************time read over******************/
        for(int i=0;i<3;i++)
        {
            data[i]=atof(str.substr(3+19+i*19,19).c_str());
        }
        for(int i=0;i<7;i++)
        {
            getline(loader,str);
            for(int j=0;j<4;j++)
            {
                // std::cout<<str.substr(3+j*19,19)<<"\n";
                data[3+j+i*4]=atof(str.substr(3+j*19,19).c_str());
            }
        }
        
        break;
    }
}

void GetBrdcPara(double *data,double *ephpara)
{
    int compasate=4;
    ephpara[0]=data[compasate+2];
    ephpara[1]=data[compasate+4];
    ephpara[2]=data[compasate+6];
    ephpara[3]=data[compasate+11];
    ephpara[4]=data[compasate+13];
    ephpara[5]=data[compasate+9];
    ephpara[6]=data[compasate+1];
    ephpara[7]=data[compasate+14];
    ephpara[8]=data[compasate+15];
    ephpara[9]=data[compasate+3];
    ephpara[10]=data[compasate+5];
    ephpara[11]=data[compasate+8];
    ephpara[12]=data[compasate+10];
    ephpara[13]=data[compasate+12];
    ephpara[14]=data[compasate+0];
    ephpara[15]=data[compasate+7];
}
