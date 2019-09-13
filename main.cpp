//implement 3D of Elasticity model
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include "io.h"
#include "mass_matrix.h"

#define MAXITER 200
#define EPSILON 1e-6
#define DELTA 0.99
#define FRAMENUM 100
#define HARDCONST 1e+100
using namespace std;
using namespace Eigen;

//global variables
static double h;
int elements, n, fix_points;
MatrixXd Dme;
MatrixXd M;//= MatrixXd::Zero();
MatrixXd Ke, K;
MatrixXd x0, x;
MatrixXd f0, fext;
MatrixXd v; //need to be thought again
MatrixXd freedom, S;

Matrix<double, -1, -1> nods;
MatrixXi tet;
Matrix<double, -1, -1> result_nod;

void mass_Be_stiff(){
    Matrix4d A, A_in;
    MatrixXd Be = MatrixXd::Zero(6,12);

    /*rearrange x0*/
    for(int j=0; j<n; j++){
        for(int t=0; t<3; t++){
            x0(3*j+t,0) = nods(t, j);
        }
        /*translation of warped model*/
        x0(3*j, 0) += 6;
    }

    //trying library
    Matrix<double, -1, 1> mass_vector;
    Matrix<double, 3, -1> nod;
    Matrix<int, 4, -1> cells;
    cells = tet;
    nod = nods;
    double rho = 1;
    marvel::calc_mass_vector(nod, cells, rho, mass_vector);
    for(int i=0; i<n; i++){
        for(int j=0; j<3; j++)
            M(3*i+j, 3*i+j) = mass_vector(i,0);
    }

    for(int i=0; i<elements; i++){
        /*compute mass*/
        //1. put vectors in matrix a
        for(int j=0; j<4; j++){
            for(int b=0; b<4; b++){
                if(b==0){
                    A(b,j) = 1;
                    continue;
                }
                A(b,j) = nods(b-1,tet(j,i));

            }
        }
        //2. count the mass in each node of an element
        double vole = fabs(A.determinant())/6;
        /*for(int j=0; j<4; j++){
            M(3*tet(j,i)  ,3*tet(j,i)  ) += (vole/4); //maybe bug :vole/12?
            M(3*tet(j,i)+1,3*tet(j,i)+1) += (vole/4);
            M(3*tet(j,i)+2,3*tet(j,i)+2) += (vole/4);
        }*/

        /*compute strain coefficient Be of e = Be*u */
        A_in = A.inverse();

        //A_in = A_in*6*vole;  //maybe bug is right here?
        for(int t=0; t<4; t++){
            for(int p=0; p<3; p++){
                //brutal force
                if(p==0) Be.col(3*t+p) << A_in(t,1),0,0,A_in(t,2),0,A_in(t,3);
                else if(p==1) Be.col(3*t+p) << 0,A_in(t,2),0,A_in(t,1),A_in(t,3),0;
                else Be.col(3*t+p) << 0,0,A_in(t,3),0,A_in(t,2),A_in(t,1);
            }
        }
        //Be = Be/(6*vole);  //to be considered;

        /*compute stiffness*/
        double young =800;
        double poisson =0.45;
        MatrixXd E = MatrixXd::Zero(6,6);{ //maybe i should use a txt file to read it in?
            E.block(0,0,3,3) << 1-poisson, poisson, poisson,
                                poisson, 1-poisson, poisson,
                                poisson, poisson, 1-poisson;
            E(3,3) = 0.5-poisson;
            E(4,4) = 0.5-poisson;
            E(5,5) = 0.5-poisson;
            E = E*young/((1+poisson)*(1-poisson*2));
        }

        MatrixXd stiff(12,12);
        stiff = vole*Be.transpose()*E*Be;

        for(int j=0; j<12; j++){
            for(int t=0; t<12; t++){
                Ke(12*i+j,12*i+t) = stiff(j,t);
            }
        }


    }
    //M = M*freedom;
    //cout << "mass:\n"<<M << endl;
}

/*class energy{
public:
    energy(const double h, MatrixXd v, MatrixXd &x)
        :h_(h), vn_(v),xn_(x){}
    void Grad(MatrixXd &grad, MatrixXd &x, MatrixXd &K )const{
        //grad += (M+h_*h_*K)*((x-xn_)/h_)+(h_*K*xn_);
        grad += (M+h_*h_*(K+S.transpose()*S*(x-x0)))*((x-xn_)/h_);
    }
    void Hess(MatrixXd &hess, MatrixXd &x, MatrixXd &K )const{
        //hess += (M+h_*h_*K)/h_;
    }
    void update(MatrixXd &x){
        //v = vn_;
        xn_ = x;
    }
private:
    const double h_;
    MatrixXd vn_, xn_;
};

class impl_euler_energy{
public:
    impl_euler_energy(MatrixXd &x,const double h, MatrixXd v)
        :xn_(x),h_(h),vn_(v){}
    void Grad(MatrixXd &grad, MatrixXd &x, MatrixXd &K )const{
        grad += -(M*vn_-h_*((K+S.transpose()*S*(x-x0))*xn_-f0*x0));
        //grad += -(M*vn_-h*K*xn_);//(M+h_*h_*K)*(x-xn_)/h_;
    }
    void Hess(MatrixXd &hess, MatrixXd &x, MatrixXd &K )const{
        hess+=S.transpose()*S;
        //hess += h_*K;
        //hess += M / (h_*h_);
    }
    void update(MatrixXd &x){
        vn_ = (x-xn_)/h_;
        v = vn_;
        xn_ = x;
    }
private:
    const double h_;
    MatrixXd vn_,xn_;
};*/
class energy{
public:
    energy(const double h, MatrixXd v, MatrixXd &x)
        :h_(h), vn_(v),xn_(x){}
    void Grad(MatrixXd &grad, MatrixXd &x, MatrixXd &K )const{
        grad += (HARDCONST*S.transpose()*S)*(x-x0)+K*x-f0*x0-fext; //co-rotation
        //grad += (K+HARDCONST*S.transpose()*S)*(x-x0)-fext; //linear
    }
    void Hess(MatrixXd &hess, MatrixXd &x, MatrixXd &K )const{
        hess+=(K+HARDCONST*S.transpose()*S);
    }
    void update(MatrixXd &x){
        xn_ = x;
    }
private:
    const double h_;
    MatrixXd vn_, xn_;
};
class impl_euler_energy{
public:
    impl_euler_energy(MatrixXd &x,const double h, MatrixXd v)
        :xn_(x),h_(h),vn_(v){}
    void Grad(MatrixXd &grad, MatrixXd &x, MatrixXd &K )const{
        grad += M*((x-xn_)/h_ - vn_)/h_;
    }
    void Hess(MatrixXd &hess, MatrixXd &x, MatrixXd &K )const{
        hess+=M/(h_*h_);
    }
    void update(MatrixXd &x){
        vn_ = (x-xn_)/h_;
        v = vn_;
        xn_ = x;
    }
private:
    const double h_;
    MatrixXd vn_,xn_;
};


impl_euler_energy *Ev;
energy *Et;

void newton(){
    MatrixXd cur_x(3*n,1);
    MatrixXd prev_x(3*n,1);
    cur_x = x;
    prev_x = cur_x;

    /*compute rotation Re*/
    MatrixXd Re = MatrixXd::Zero(12*elements, 12*elements);
    //compute displacement gradient
    MatrixXd Dse = MatrixXd::Zero(3*elements,3*elements);

    for(int i=0; i<elements; i++){
        for(int j=0; j<3; j++){
            for(int b=0; b<3; b++){
                Dse(3*i+b,3*i+j) = cur_x(3*tet(j,i)+b,0)-cur_x(3*tet(3,i)+b,0);
            }
        }
    }

    for(int i=0; i<elements; i++){
        //get Dm and Ds
        Matrix3d Ds, Dm;
        for(int j=0; j<3; j++){
            for(int t=0; t<3; t++){
                Ds(t,j) = Dse(3*i+t, 3*i+j);
                Dm(t,j) = Dme(3*i+t, 3*i+j);
            }
        }

        MatrixXd rot = MatrixXd::Zero(3, 3);
        rot = Ds*Dm.inverse(); //initial rot_matrix with displacement gradient
        //if(i==0) cout <<"rot: "<<rot<<endl;
        JacobiSVD<MatrixXd> svd(rot, ComputeFullU | ComputeFullV);
        MatrixXd rotation = MatrixXd::Zero(3,3);
        rotation = svd.matrixU()*svd.matrixV().transpose();


        for(int j=0; j<3; j++){
            for(int t=0; t<3; t++){
                Re(12*i+j,12*i+t) = rotation(j,t);
                Re(12*i+j+3,12*i+t+3) = rotation(j,t);
                Re(12*i+j+6,12*i+t+6) = rotation(j,t);
                Re(12*i+j+9,12*i+t+9) = rotation(j,t);
            }
        }
    }

    /*stiffness*/
    MatrixXd R(12,12), stiff(12,12);
    K = MatrixXd::Zero(3*n, 3*n);
    f0 = MatrixXd::Zero(3*n, 3*n);
    for(int i=0; i<elements; i++){
        for(int j=0; j<12; j++){
            for(int t=0; t<12; t++){
                R(j,t) = Re(12*i+j,12*i+t);
                stiff(j,t) = Ke(12*i+j,12*i+t);
            }
        }

        MatrixXd tmp_k = R*stiff*R.transpose();

        for(int j=0; j<4; j++){
            for(int t=0; t<4; t++){
                for(int p=0; p<3; p++){
                    for(int q=0; q<3; q++){
                        K(3*tet(j,i)+p,3*tet(t,i)+q) += tmp_k(3*j+p, 3*t+q);
                    }
                }
            }
        }
    /*compute initial force*/
        MatrixXd tmp_f;
        tmp_f = R*stiff;//*tmp_x;

        for(int j=0; j<4; j++){
            for(int t=0; t<4; t++){
                for(int p=0; p<3; p++){
                    for(int q=0; q<3; q++){
                        f0(3*tet(j,i)+p,3*tet(t,i)+q) += tmp_f(3*j+p, 3*t+q);//think again
                        //f0(3*tet(j,i)+t,0) += tmp_f(j*3+t,0);
                    }
                }
            }
        }
    }
    //f0 = K*x0;
    //f0 = f0*(-1);

    /*external force*/
    for(int i=0; i<n; i++){
        fext(3*i,0) = -9.8*5; //gravity
        fext = M*fext;
    }

    /*solving equation*/
    for(int i=0; i<MAXITER; i++){
        MatrixXd grad = MatrixXd::Zero(3*n, 1);
        Ev->Grad(grad, cur_x, K);
        Et->Grad(grad, cur_x, K);
        grad *= (-1);
        cout << "grad: " <<grad.norm() << endl;

        MatrixXd hess = MatrixXd::Zero(3*n, 3*n);
        Ev->Hess(hess, cur_x, K);
        Et->Hess(hess, cur_x, K);
        cout << "hess: " <<hess.norm() << endl;

        //determine convergence condition
        MatrixXd dx = hess.inverse()*grad;

        cout <<"current: " << cur_x.norm() << " dx: " <<dx.norm() <<" prev:"<< prev_x.norm() << endl;
        if(dx.norm()<= EPSILON*prev_x.norm()){
            cout << "\t@CONVERGED\n";
            break;
        }
        //update current position
        prev_x = cur_x;
        cur_x = cur_x + dx;

    }
    x = cur_x;
}

void transfer(){
    result_nod = nods;
    for(int i=0; i<n; i++){
        for(int j=0; j<3; j++){
            result_nod(j,i) = x(3*i+j,0);
        }
    }
}

int main(){
    char str[200] = "C:\\Users\\iivy303\\Desktop\\experi\\frames\\frm_%d.vtk";
    //input: extracts the matrix from .vtk
    marvel::tet_mesh_read_from_vtk("beam-0.1k.vtk", nods, tet);
    n = nods.cols();    //suppose not a scalar
    fix_points = 0;

    /*pre-computation Be, Ke, M, decomposition gradient's denominator*/
    elements = tet.cols();
    x0.resize(3*n,1);
    x.resize(3*n,1);
    v = MatrixXd::Zero(3*n,1); //carefully think it again
    //f0 = MatrixXd::Zero(3*n, 1);
    f0 = MatrixXd::Zero(3*n, 3*n); //exp
    fext = MatrixXd::Zero(3*n,1);
    Dme = MatrixXd::Zero(3*elements,3*elements);
    K = MatrixXd::Zero(3*n, 3*n);
    Ke = MatrixXd::Zero(12*elements,12*elements);
    M = MatrixXd::Zero(3*n, 3*n);
    freedom = MatrixXd::Identity(3*n,3*n);


    for(int i=0; i<elements; ++i){
        for(int j=0; j<3; j++){
            for(int b=0; b<3; b++){
                Dme(3*i+b,3*i+j) = nods(b,tet(j,i))-nods(b,tet(3,i));
                //cout << "dme:\n" << Dme(3*i+b,3*i+(j-1)) << endl;
            }
        }
    }

    mass_Be_stiff();

    /*initialize*/
    x = x0; //position
    h = 0.1; //time step
    /*//circumstance#1: giving some initial velocity
    for(int i=0; i<n; i++){
        v(3*i+1,0) = pow(x0(3*i,0),2)-12; //velocity of y-direction = x^2 -12
        if(x0(3*i+1,0)<0.9){
            //v(3*i+1,0) = 0;
            freedom(3*i+1,3*i+1) = 1.7e+200;
        }
    }*/
    //circumstance#2: fix one end of the rod by selection matrix
    MatrixXd tmp_fix = nods.row(2);
    double mini = tmp_fix.minCoeff();
    cout <<"mini: "<<mini <<endl;
    for(int i=0; i<n; i++){
        if(x0(3*i+2,0) < mini*DELTA){
            fix_points++;
        }
    }
    S = MatrixXd::Zero(3*fix_points, 3*n);
    int j=0;
    for(int i=0; i<n; i++){
        if(x0(3*i+2,0) < mini*DELTA){ //think again
            S(3*j, 3*i) = 1;
            S(3*j+1, 3*i+1) = 1;
            S(3*j+2, 3*i+2) = 1;
            j++;
        }
        if(j>=fix_points) break;
    }


    Ev = new impl_euler_energy(x, h, v);
    Et = new energy(h, v, x);

    double ek[FRAMENUM], u[FRAMENUM], tot[FRAMENUM];
    //simulate
    for(int frame=0; frame<FRAMENUM; frame++){
        cout << "frm: "<< frame << " position: "<<x.norm()<<endl;
        newton();

        Ev->update(x);
        Et->update(x);

        /*check energy if it is stable*/
        MatrixXd e_kinectic = 0.5*v.transpose()*M*v;
        MatrixXd tmp = (x-x0);
        MatrixXd u_energy = 0.5*tmp.transpose()*K*tmp; //suppose .cwiseAbs() ?

        ek[frame] = e_kinectic(0,0);
        u[frame] = u_energy(0,0);
        tot[frame] = e_kinectic(0,0) + u_energy(0,0);

        /*output: output to file*/
        transfer();
        //sprintf(str,"C:\\Users\\iivy303\\Desktop\\experi\\frames\\frm_%d.vtk",frame);
        sprintf(str,"C:\\Users\\iivy303\\Desktop\\experi\\temp\\rot_%d.vtk",frame);
        //ofstream out(str);
        marvel::tet_mesh_write_to_vtk<double>(str, result_nod, tet);
        //out.close();
    }
    cout << "[DONE]\n";

    fstream file;
    file.open("energy.txt",ios::out);
    for(int i=0; i<FRAMENUM; i++){
        file<<i    <<"  ";
        file<<ek[i]<<"   ";
        file<<u[i]<<"   ";
        file<<tot[i]<<"\n";

    }
    file.close();

    return 0;
}

