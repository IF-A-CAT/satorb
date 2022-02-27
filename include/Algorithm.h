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
    private:
    Matrix Pxx_post;                                //验后
    double corr;                                    //单位权中误差
    Matrix v;                                       //残差
    Matrix x;                                        //parameter
    Matrix Pxx_info;                                //先验
    Matrix l;                                   
    Matrix B;                                        //designed matrix
};

void Split(std::string str,const char sign,double * data); 


#endif