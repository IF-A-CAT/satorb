#pragma once

#ifndef _ALGORITHM_H
#define _ALGORITHM_H
#include"VecMat.h"
//  This demo is for algorithm used in GNSS  //
class Lsq
{   public:
    Lsq(const Matrix &Pxx_info,const Matrix &l,const Matrix &B,int n);
    Lsq(const Matrix &l,const Matrix &B,int n);
    Matrix get_pxx(){return Pxx_post;};
    double get_corr(){return corr;};
    Matrix get_x(){return x;};
    Matrix get_v(){return v;};
    private:
    Matrix Pxx_post;                                //验后
    double corr;                                    //单位权中误差
    Matrix v;                                       //残差
    Matrix x;                                        //parameter
    Matrix Pxx_info;                                //先验
    Matrix l;                                   
    Matrix B;                                        //designed matrix
};


class EKF
{   public:
    EKF(int n_);

    void init(const Matrix &x,const Matrix &P);
    void init(int epoch_,double time_ ,const Matrix &x,const Matrix &P);

    void TimeUpdate(double time,const Matrix & X_,const Matrix &STM);               //STM: State Tranform Matrix

    void MeasUpdate(const Matrix &Z,const Matrix &g,const Matrix Weight_inv,const Matrix &G);

    Matrix get_X()const{return X;};
    Matrix get_Pxx()const{return Pxx;};
    double get_time()const {return time;};
    int  get_epoch()const{return epoch;};

    private:
    Matrix X;     //state vector
    Matrix Pxx;   //number of state parameter
    int epoch ;    
    int n;
    double time;
    
};


void get_coef_mat(const Matrix &obs,const Matrix &ref,double clk,Matrix &H,Matrix &L);

bool gross_error_detection(int freedom,double sigma,const Matrix &v,const Matrix &Weight);

void Split(std::string str,const char sign,double * data); 




#endif