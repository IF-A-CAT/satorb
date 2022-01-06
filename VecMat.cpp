#include"VecMat.h"
#include<iostream>
#include<cmath>
#include<iomanip>
using namespace std;

Vector::Vector():n(0){};

Vector::Vector(int i):n(i){
    V=new double [n];
};

Vector::Vector(double *v,int n):n(n){
    V=new double [n];
    for(int i=0;i<n;i++)    V[i]=v[i];
}

Vector::Vector(double a,double b,double c,double d,double e,double f):n(6){
    V=new double [6];
    V[0]=a;V[1]=b;V[2]=c;V[3]=d;V[4]=e;V[5]=f;
}

Vector::Vector(double a,double b,double c):n(3){
    V=new double [3];
    V[0]=a;V[1]=b;V[2]=c;
}

Vector::Vector(const Vector &V1):n(V1.n){
    V=new double [n];
    for(int i=0;i<n;i++)    V[i]=V1.V[i];
}


Vector::~Vector(){
    delete [] V;
}


Vector& Vector::operator=(const double &value){
    n=1;
    V=new double [1];
    V[0]=value;
    return *this;
}

Vector& Vector::operator=(const Vector &V1){
    n=V1.n;
    V=new double [n];
    for(int i=0;i<n;i++)    V[i]=V1.V[i];
    return *this;
}

Vector operator+(const Vector &V1,const Vector &V2){
    if(V1.n!=V2.n){
        cerr<<"ERROR:Vector 1 and Vector 2's dims are different,please check them"<<endl;
        exit(1);}
    Vector V3(V1.n);
    for(int i=0;i<V1.n;i++) V3.V[i]=V2.V[i]+V1.V[i];
    return V3;
}

Vector operator+(const double &value,const Vector &V){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=value+V.V[i];
    return V1;
}

Vector operator+(const Vector &V,const double &value){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=value+V.V[i];
    return V1;
}

Vector operator-(const Vector &V1,const Vector &V2){
    if(V1.n!=V2.n){
        cerr<<"ERROR:Vector 1 and Vector 2's dims are different,please check them"<<endl;
        exit(1);}
    Vector V3(V1.n);
    for(int i=0;i<V1.n;i++) V3.V[i]=-V2.V[i]+V1.V[i];
    return V3;
}

Vector operator-(const double &value,const Vector &V){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=value-V.V[i];
    return V1;
}

Vector operator-(const Vector &V,const double &value){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=-value+V.V[i];
    return V1;
}

Vector operator-(const Vector &V){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=-V.V[i];
    return V1;
}

Vector operator/(const Vector &V,const double &value){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=V.V[i]/value;
    return V1;
}

Vector operator*(const Vector &V,const double &value){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=V.V[i]*value;
    return V1;
}

Vector operator*(const double &value,const Vector &V){
    Vector V1(V.n);
    for(int i=0;i<V1.n;i++) V1.V[i]=V.V[i]*value;
    return V1;
}

double Norm(const Vector &V){
    double length=0;
    for(int i=0;i<V.n;i++)    length+=pow(V.V[i],2);
    return sqrt(length);
}

double Dot(const Vector &V1,const Vector &V2){
    if(V1.n!=V2.n){
        cerr<<"ERROR:Vector 1 and Vector 2's dims are different,please check them"<<endl;
        exit(1);}
    double dot_result=0;
    for(int i=0;i<V1.n;i++) dot_result+=V1.V[i]*V2.V[i];
    return dot_result;
}

Vector Cross3D(const Vector &V1,const Vector &V2){
    if(V1.n!=3||V2.n!=3){
        cerr<<"ERROR:Vectors must be 3-dims,please check them"<<endl;
        exit(1);}
    Vector V3(3);
    V3.V[0]=V1.V[1]*V2.V[2]-V1.V[2]*V2.V[1];
    V3.V[1]=V1.V[2]*V2.V[0]-V1.V[0]*V2.V[2];
    V3.V[2]=V1.V[0]*V2.V[1]-V1.V[1]*V2.V[0];
    return V3;
}

ostream& operator<<(ostream &os,const Vector &V){
    int width=os.width();
    for(int i=0;i<V.n;i++)
        os<<setw(8)<<V(i)<<endl;
    return os;
}

Matrix::Matrix():r(0),c(0){
    M=new double *[1];
}

Matrix::Matrix(int R,int C):r(R),c(C){
    M=new double*[r];
    int i,j;
    for(i=0;i<r;i++)    M[i]=new double [c];
    for(i=0;i<r;i++){
        for(j=0;j<c;j++)
            M[i][j]=0;
    }
}
Matrix::Matrix(const Matrix &M1):c(M1.c),r(M1.r){
    M=new double*[r];
    int i,j;
    for(i=0;i<r;i++)    M[i]=new double [c];
    for(i=0;i<r;i++){
        for(j=0;j<c;j++)
            M[i][j]=M1.M[i][j];
    }
}

Matrix::Matrix(const Vector &V):r(V.n),c(1){
    M=new double *[r];
    int i;
    for(i=0;i<r;i++){
        M[i]=new double [1];
        M[i][0]=V.V[i];
    }
}

Matrix::Matrix(double *m,int R,int C):r(R),c(C){
    M=new double*[r];
    int i;
    for(i=0;i<r;i++)    M[i]=new double [c];
    for(i=0;i<r;i++){
        for(int j=0;j<c;j++)
            M[i][j]=m[i*c+j];
    }
}

Matrix& Matrix::operator=(const double &value){
    M=new double *[1];
    M[0]=new double [1];
    M[0][0]=value;
    r=1;
    c=1;
    return *this;
}

Matrix::~Matrix(){
    delete [] M;
}

Matrix& Matrix::operator=(const Vector &V){
    M=new double*[V.n];
    r=V.n;c=1;
    int i;
    for(i=0;i<r;i++)    {
        M[i]=new double [1];
        M[i][0]=V.V[i];
    }
    return *this;
}

Matrix operator+(const Matrix & M,const double &value){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=M(i,j)+value;
    }
    return M1;
}
Matrix operator+(const double &value,const Matrix & M){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=M(i,j)+value;
    }
    return M1;
}

Matrix operator+(const Matrix & M1,const Matrix &M2){
    if(M1.r!=M2.r||M1.c!=M2.c){
        cerr<<"ERROR:Tow matrixes' dims are different,please check them"<<endl;
        exit(1);
    }
    Matrix M3(M1);
    for(int i=0;i<M3.r;i++){
        for(int j=0;j<M3.c;j++)
            M3(i,j)=M2(i,j)+M1(i,j);
    }
    return M3;
}

Matrix operator-(const Matrix & M,const double &value){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M.c;j++)
            M1(i,j)=M(i,j)-value;
    }
    return M1;
}

Matrix operator-(const double &value,const Matrix & M){
   Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M.c;j++)
            M1(i,j)=-M(i,j)+value;
    }
    return M1; 
}
Matrix operator-(const Matrix & M1,const Matrix &M2){
    if(M1.r!=M2.r||M1.c!=M2.c){
        cerr<<"ERROR:Tow matrixes' dims are different,please check them"<<endl;
        exit(1);
    }
    Matrix M3(M1);
    for(int i=0;i<M3.r;i++){
        for(int j=0;j<M3.c;j++)
            M3(i,j)=M1(i,j)-M2(i,j);
    }
    return M3;
}
Matrix operator-(const Matrix &M){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=-M(i,j);
    }
    return M1;
}

Matrix operator*(const double &value,const Matrix &M){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=value*M(i,j);
    }
    return M1;
}
Matrix operator*(const Matrix & M,const double &value){
    Matrix M1(M);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=value*M(i,j);
    }
    return M1;
}
Matrix operator*(const Matrix & M1,const Matrix &M2){
    if(M1.c!=M2.r){
        cerr<<"ERROR:The column of Matrix 1 is different from the row of Matrix 2"<<endl;
        exit(1);
    }
    Matrix M3(M1.r,M2.c);
    int i,j,k;
    for(i=0;i<M3.r;i++){
        for(j=0;j<M3.c;j++)
            for(k=0;k<M1.c;k++)
                M3(i,j)+=M1(i,k)*M2(k,j);
    }
    return M3;
}
Matrix operator*(const Matrix & M,const Vector &V){
    if(M.c!=V.n){
        cerr<<"ERROR:The colum of the matrix is different from the dim of vector"<<endl;
        exit(1);
    }
    Matrix M1(M.r,1);
    for(int i=0;i<M.r;i++){
        for(int j=0;j<V.n;j++)
            M1(i,0)+=M(i,j)*V(j);
    }
    return M1;
}

Matrix operator*(const Vector & V,const Matrix &M){
    if(M.r!=1){
        cerr<<"ERROR:The row of the matrix is not 1"<<endl;
        exit(0);
    }
    Matrix M1(V.n,M.c);
    for(int i=0;i<V.n;i++){
        for(int j=0;j<M.c;j++)
            M1(i,j)=V(i)*M(0,j);
    }
    return M1;
}


Matrix operator/(const Matrix &M,const double &value){
    Matrix M1(M);
    for (int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=M(i,j)/value;
    }
    return M1;
}

Matrix& Matrix::operator=(const Matrix &M1){
    r=M1.r;
    c=M1.c;
    M=new double*[r];
    int i,j;
    for(i=0;i<r;i++)    M[i]=new double [c];
    for(i=0;i<r;i++){
        for(j=0;j<c;j++)
            M[i][j]=M1.M[i][j];
    }
    return *this;
}

ostream& operator<<(ostream &os,const Matrix &M){
    for(int i=0;i<M.r;i++){
        for(int j=0;j<M.c;j++)
            if(j!=M.c-1){
                os<<setw(12)<<M(i,j)<<',';
            }
            else{
                os<<setw(12)<<M(i,j)<<endl;
            }
    }
    return os;
}

Matrix Eye(int n){
    Matrix eye(n,n);
    for(int i=0;i<n;i++){
        eye(i,i)=1;
    }
    return eye;
}

Matrix Trans(const Matrix &M){
    Matrix M1(M.c,M.r);
    for(int i=0;i<M1.r;i++){
        for(int j=0;j<M1.c;j++)
            M1(i,j)=M(j,i);
    }
    return M1;
}

Matrix Matrix::slice(int r1,int r2,int c1,int c2)const {
    if(r1<0||r2<0||c1<0||c2<0||r2>r-1||c2>c-1||r1>r2||c1>c2){
        cerr<<"ERROR:Index is out of range!"<<endl;
        exit(1);
    }
    int R=r2-r1+1,C=c2-c1+1;
    Matrix back(R,C);
    for(int i=0;i<R;i++){
        for(int j=0;j<C;j++)
            back(i,j)=M[r1+i][c1+j];
    }
    return back;
}

double Norm(const Matrix &M){
    if(M.c!=1){
        cerr<<"ERROR:Matrix's column is not 1,Norm fial"<<endl;
        exit(1);
    }
    double norm=0;
    for(int i=0;i<M.r;i++){
        norm+=M(i,0)*M(i,0);
    }
    return sqrt(norm);
}

// Matrix Inv_LU(const Matrix &M){

// }
Matrix* QR_decompose(const Matrix &M){
    int r=M.r,c=M.c;
    double norm;
    Matrix w(r,1),Aux(r,1),trans(r,r),M1(M);
    Matrix* QR=new Matrix [2];
    QR[0]=Eye(r);
    for(int i=0;i<r&&i<c;i++){
        Matrix  e(r-i,1);
        for(int j=0;j<r-i;j++){
            if(j!=0){
                e(j,0)=0;
            }
            else{
                e(j,0)=1;
            }
        }
        Aux=M1.slice(i,r-1,i,i);
        norm=Norm(Aux);
        if(norm<eps){
            cerr<<"ERROR:The matrix is not column full rank,or the square matrix is singular"<<endl;
            exit(1);
        }
        // cout<<Aux;
        Aux=Aux/norm;
        // cout<<Aux;
        Aux=Aux-e;
        // cout<<Aux;
        norm=Norm(Aux);
        if(norm<eps){
            continue;
        }
        Aux=Aux/norm;
        trans=Eye(r-i);

        Aux=(trans-2*Aux*Trans(Aux));
        // cout<<Aux*Trans(Aux);
        // std::cout<<Aux<<Trans(Aux)<<endl;
        trans=Eye(r);
        for(int k=0;k<r-i;k++){
            for(int m=0;m<r-i;m++)
                trans(i+k,m+i)=Aux(k,m);
        }
        QR[0]=trans*QR[0];
        M1=trans*M1;
    }
    QR[0]=Trans(QR[0]);
    QR[1]=M1;
    return QR;
   
}
// Matrix Inv_QR(const Matrix &M){
    // if(M.c!=M.r){
        // cerr<<"ERROR:Matrix is not a square matrix"<<endl;
        // exit(1);
    // }
 
// }