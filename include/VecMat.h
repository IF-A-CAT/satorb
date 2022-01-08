#pragma once 

#ifndef _VECMAT_H
#define _VECMAT_H
//This for matrix calculation

// #include<iostream>
#include<iostream>

static double eps=__DBL_EPSILON__*100;

class Matrix;


class Vector
{
    public:
     friend class Matrix;
    //init
     Vector();
     Vector(const Vector &V );             //create vectorusing vector
     Vector(double a,double b,double c);   // 3-dim vector
     Vector(double a,double b,double c,double d,double e,double f);//6-dim
     Vector(double *v,int n);                  //set vector with array
     Vector(int n);                      //dim?
     
 
     ~Vector();
     //compute
     double operator()(int n) const{return V[n];};
     double& operator()(int n){return V[n];};
     
     Vector& operator=(const Vector &V);
     Vector& operator=(const double &v);
     friend Vector operator+(const Vector &V1,const Vector &V2);
     friend Vector operator+(const double &v,const Vector &V);
     friend Vector operator+(const Vector &V,const double &v);
 
 
     friend Vector operator-(const Vector &V1,const Vector &V2);
     friend Vector operator-(const double &v,const Vector &V);
     friend Vector operator-(const Vector &V,const double &v);
     friend Vector operator-(const Vector &V);

     friend Vector operator/(const Vector &V,const double &v);
 
     friend Vector operator*(const Vector &V,const double &v);
     friend Vector operator*(const double &v,const Vector &V);
     
     friend Matrix operator*(const Matrix & M,const Vector &V);
     friend Matrix operator*(const Vector & V,const Matrix &M);

     friend double Norm(const Vector &V);
     friend double Dot(const Vector &V1,const Vector &V2);
     friend Vector Cross3D(const Vector &V1,const Vector &V2);
     friend std::ostream& operator<<(std::ostream &os,const Vector &V);

    private:
     double *V;    //vector pointer
     int n;    //The length of the vector
};



class Matrix
{
    public:
        Matrix();
        Matrix(int r,int c);
        Matrix(double *array,int r,int c);
        Matrix(const Matrix &M);
        Matrix(const Vector &V);

        ~Matrix();

        double operator()(int r,int c) const{return M[r][c];};
        double& operator()(int r,int c) {return M[r][c];};

        Matrix& operator=(const double &value);
        Matrix& operator=(const Vector &V);
        Matrix& operator=(const Matrix &M);

        friend Matrix operator+(const Matrix & M,const double &value);
        friend Matrix operator+(const double &value,const Matrix & M);
        friend Matrix operator+(const Matrix & M1,const Matrix &M2);

        friend Matrix operator-(const Matrix & M,const double &value);
        friend Matrix operator-(const double &value,const Matrix & M);
        friend Matrix operator-(const Matrix & M1,const Matrix &M2);
        friend Matrix operator-(const Matrix &M);

        friend Matrix operator*(const double &value,const Matrix &M);
        friend Matrix operator*(const Matrix & M,const double &value);
        friend Matrix operator*(const Matrix & M1,const Matrix &M2);
        friend Matrix operator*(const Matrix & M,const Vector &V);
        friend Matrix operator*(const Vector & V,const Matrix &M);

        friend Matrix* QR_decompose(const Matrix &M);
        friend Matrix operator/(const Matrix &M,const double &value);








        
        friend Matrix Eye(int n);
        friend Matrix Trans(const Matrix &M);
        friend Matrix Inv_LU(const Matrix &M); 
        friend Matrix Inv_QR(const Matrix &M);
        friend Matrix LUinv(const Matrix &M);
        friend double Norm(const Matrix &M);
        friend bool Is_L(const Matrix &M);
        Matrix slice(int r1,int r2,int c1,int c2)const;
        friend std::ostream& operator<<(std::ostream &os,const Matrix &V);
        friend Matrix R_x( double alpha);
        friend Matrix R_y( double alpha);
        friend Matrix R_z( double alpha);

    private:
        double **M;
        int r;
        int c;
};
#endif