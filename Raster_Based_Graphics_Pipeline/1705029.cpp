#include<bits/stdc++.h>
#include<iostream>
#include<fstream>
#include<string>
#include <sstream>
#include<stack>
#include<vector>
#include<iomanip>
#include <cstdlib>
#define pi 2.0*acos(0.0)
#define INF numeric_limits<double>::infinity()
#include "bitmap_image.hpp"
using namespace std;
double eyeX,eyeY,eyeZ;
double lookX,lookY,lookZ;
double upX,upY,upZ;
double fovY,aspectRatio,near,far;
int screenWidth,screenHeight;
double left_lim_X,right_lim_X,bottom_lim_Y,top_lim_Y,frontZ,rearZ;
ofstream st_1;
ofstream st_2;
ofstream st_3;
ofstream st_4;
struct point
{
    double x;
    double y;
    double z;
    double w;
};
struct color
{
    int r;
    int g;
    int b;

};
void printMatrix(double matrix[4][4])
{
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            cout<<matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}
double getMax(double a, double b,double c)
{
    return(a>b?(a>c?a:c):(b>c?b:c));
}
double getMin(double a, double b,double c)
{
    return(a<b?(a<c?a:c):(b<c?b:c));
}
double getintersection(double x1,double y1,double x2,double y2,double scan_line)
{
    double t=(scan_line-y1)/(y2-y1);
    double multiplied=t*(x2-x1);
    return (x1+multiplied);

}
bool isInvalid(point p,double x1,double x2,double y1,double y2)
{
    if(p.x>max(x1,x2)||p.x<min(x1,x2)||p.y>max(y1,y2)|p.y<min(y1,y2)) return true;
    else return false;
}

point i,j,k,look,eye,up;
double dot(point a,point b)
{
    double p=(a.x*b.x)+(a.y*b.y)+(a.z*b.z);
    return p;
}



point cross(point p1,point p2)
{
    point result;
    result.x=(p1.y*p2.z)-(p2.y*p1.z);
    result.y=(p1.z*p2.x)-(p1.x*p2.z);
    result.z=(p1.x*p2.y)-(p1.y*p2.x);
    return result;
}
point normalize(point p)
{
    point temp;
    double r=sqrt((p.x*p.x)+(p.y*p.y)+(p.z*p.z));
    temp.x=p.x/r;
    temp.y=p.y/r;
    temp.z=p.z/r;
    return temp;
}
point sub(point p1,point p2)
{
    point result;
    result.x=p1.x-p2.x;
    result.y=p1.y-p2.y;
    result.z=p1.z-p2.z;
    return result;
}

point Rodrigues(point axis,point a,double angle)
{
    point rodri;
    double d=dot(a,axis);
    point c=cross(a,axis);
    rodri.x=(cos((angle*pi)/180.0)*axis.x)+((1-cos((angle*pi)/180.0))*d*a.x)+(sin((angle*pi)/180.0)*c.x);
    rodri.y=(cos((angle*pi)/180.0)*axis.y)+((1-cos((angle*pi)/180.0))*d*a.y)+(sin((angle*pi)/180.0)*c.y);
    rodri.z=(cos((angle*pi)/180.0)*axis.z)+((1-cos((angle*pi)/180.0))*d*a.z)+(sin((angle*pi)/180.0)*c.z);
    return rodri;

}
class transMatrix
{

public:
    double matrix[4][4];
    transMatrix()
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                if(i==j) matrix[i][j]=1;
                else matrix[i][j]=0;
            }
        }
    }
    void setMatrix(double mat[4][4])
    {
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; i++)
            {
                matrix[i][j]=mat[i][j];
            }
        }
    }
    void generate_scaling_matrix(point p)
    {
        matrix[0][0]=p.x;
        matrix[1][1]=p.y;
        matrix[2][2]=p.z;

    }
    void generate_translation_matrix(point p)
    {
        matrix[0][3]=p.x;
        matrix[1][3]=p.y;
        matrix[2][3]=p.z;

    }
    void generate_rotation_matrix(point p,double angle)
    {
        point norm,c1,c2,c3;
        norm=normalize(p);
        c1=Rodrigues(i,norm,angle);
        c2=Rodrigues(j,norm,angle);
        c3=Rodrigues(k,norm,angle);

        matrix[0][0]=c1.x;
        matrix[0][1]=c2.x;
        matrix[0][2]=c3.x;

        matrix[1][0]=c1.y;
        matrix[1][1]=c2.y;
        matrix[1][2]=c3.y;

        matrix[2][0]=c1.z;
        matrix[2][1]=c2.z;
        matrix[2][2]=c3.z;
    }
    void generate_rotation_matrix_for_view(point l,point r,point u)
    {
        matrix[0][0]=r.x;
        matrix[0][1]=r.y;
        matrix[0][2]=r.z;

        matrix[1][0]=u.x;
        matrix[1][1]=u.y;
        matrix[1][2]=u.z;

        matrix[2][0]=-l.x;
        matrix[2][1]=-l.y;
        matrix[2][2]=-l.z;
    }
    void generate_projection_matrix(double fovX)
    {


        double r = near*tan(((fovX*0.5)*pi)/180.0);
        double f=far-near;
        double t= near*tan(((fovY*0.5)*pi)/180.0);
        matrix[0][0]=near/r;
        matrix[1][1]=near/t;
        matrix[2][2]=-(far+near)/f;
        matrix[2][3]=-(2*far*near)/f;
        matrix[3][2]=-1.0;
        matrix[3][3]=0.0;


    }
    void matrix_multiplication(double mat1[4][4],double mat2[4][4])
    {
        //double temp[4][4];
        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                matrix[i][j] = 0.0;

                for(int k=0; k<4; k++)
                {
                    matrix[i][j] += mat1[i][k]*mat2[k][j];
                }
            }
        }
        //setMatrix(temp);
    }
};

class Triangle
{
public:
    point points[3];
    color colors;
};
stack <transMatrix> matrixstack;
stack <transMatrix> push_pop_stack;


point matrix_multi(double mat[4][4],point p)
{
    point temp;


    for(int i=0; i<4; i++)
    {
        double tempi=0.0;
        for(int j=0; j<4; j++)
        {
            if(j==0) tempi=tempi+(mat[i][j]*p.x);

            else if(j==1) tempi=tempi+(mat[i][j]*p.y);
            else if(j==2) tempi=tempi+(mat[i][j]*p.z);
            else tempi=tempi+(p.w*mat[i][j]);
        }

        if(i==0) temp.x=tempi;
        else if(i==1) temp.y=tempi;
        else if (i==2) temp.z=tempi;
        else
        {

            temp.w=tempi;

        }
    }
    return temp;
}
point scale(point p)
{
    point temp;

    temp.x=p.x/p.w;
    temp.y=p.y/p.w;
    temp.z=p.z/p.w;
    return temp;
}
void print_triangle(point p1,point p2,point p3,int flag)
{


    if(flag==1)
    {
        st_1<<fixed<<setprecision(7)<<p1.x;
        st_1<<fixed<<setprecision(7)<<" "<<p1.y;
        st_1<<fixed<<setprecision(7)<<" "<<p1.z<<endl;
        st_1<<fixed<<setprecision(7)<<p2.x;
        st_1<<fixed<<setprecision(7)<<" "<<p2.y;
        st_1<<fixed<<setprecision(7)<<" "<<p2.z<<endl;
        st_1<<fixed<<setprecision(7)<<p3.x;
        st_1<<fixed<<setprecision(7)<<" "<<p3.y;
        st_1<<fixed<<setprecision(7)<<" "<<p3.z<<endl;
        st_1<<endl;
    }
    else if(flag==2)
    {
        st_2<<fixed<<setprecision(7)<<p1.x;
        st_2<<fixed<<setprecision(7)<<" "<<p1.y;
        st_2<<fixed<<setprecision(7)<<" "<<p1.z<<endl;
        st_2<<fixed<<setprecision(7)<<p2.x;
        st_2<<fixed<<setprecision(7)<<" "<<p2.y;
        st_2<<fixed<<setprecision(7)<<" "<<p2.z<<endl;
        st_2<<fixed<<setprecision(7)<<p3.x;
        st_2<<fixed<<setprecision(7)<<" "<<p3.y;
        st_2<<fixed<<setprecision(7)<<" "<<p3.z<<endl;
        st_2<<endl;
    }
    else
    {
        st_3<<fixed<<setprecision(7)<<p1.x;
        st_3<<fixed<<setprecision(7)<<" "<<p1.y;
        st_3<<fixed<<setprecision(7)<<" "<<p1.z<<endl;
        st_3<<fixed<<setprecision(7)<<p2.x;
        st_3<<fixed<<setprecision(7)<<" "<<p2.y;
        st_3<<fixed<<setprecision(7)<<" "<<p2.z<<endl;
        st_3<<fixed<<setprecision(7)<<p3.x;
        st_3<<fixed<<setprecision(7)<<" "<<p3.y;
        st_3<<fixed<<setprecision(7)<<" "<<p3.z<<endl;
        st_3<<endl;
    }

    //st_1<<p1.y<<endl;
    //cout<<p1.z<<endl;

    //st_1.close();

}
int main ()
{

    string lineOfFile;
    vector<Triangle> triangles;
    transMatrix tm;
    i.x=1;
    i.y=0;
    i.z=0;
    j.x=0;
    j.y=1;
    j.z=0;
    k.x=0;
    k.y=0;
    k.z=1;
    st_1.open("stage_1.txt");
    int stack_size=0;

    matrixstack.push(tm);

    ifstream filein("scene.txt");
    if (!filein)
    {
        cout<<"Error in opening scene.txt"<<endl;
        exit(1);
    }
    else
    {
        int linecount=1;
        while (!filein.eof())
        {
            // Output the text from the file
            if(linecount==1)
            {


                getline(filein,lineOfFile);
                linecount++;

                std::istringstream iss(lineOfFile);
                iss>>eyeX>>eyeY>>eyeZ;
                eye.x=eyeX;
                eye.y=eyeY;
                eye.z=eyeZ;
            }
            if(linecount==2)
            {
                getline(filein,lineOfFile);
                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>lookX>>lookY>>lookZ;
                look.x=lookX;
                look.y=lookY;
                look.z=lookZ;

            }
            if(linecount==3)
            {
                getline(filein,lineOfFile);
                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>upX>>upY>>upZ;
                up.x=upX;
                up.y=upY;
                up.z=upZ;

            }
            if(linecount==4)
            {
                getline(filein,lineOfFile);
                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>fovY>>aspectRatio>>near>>far;

            }
            else
            {
                getline(filein,lineOfFile);
                linecount++;
                if(lineOfFile.compare("triangle")==0)
                {

                    double point_x,point_y,point_z;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss(lineOfFile);
                    iss>>point_x>>point_y>>point_z;
                    struct point point_1;
                    point_1.x=point_x;
                    point_1.y=point_y;
                    point_1.z=point_z;
                    point_1.w=1.0;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss1(lineOfFile);
                    iss1>>point_x>>point_y>>point_z;
                    struct point point_2;
                    point_2.x=point_x;
                    point_2.y=point_y;
                    point_2.z=point_z;
                    point_2.w=1.0;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss2(lineOfFile);
                    iss2>>point_x>>point_y>>point_z;
                    struct point point_3;
                    point_3.x=point_x;
                    point_3.y=point_y;
                    point_3.z=point_z;
                    point_3.w=1.0;
                    point p1,p2,p3;

                    p1=matrix_multi(matrixstack.top().matrix,point_1);
                    p2=matrix_multi(matrixstack.top().matrix,point_2);
                    p3=matrix_multi(matrixstack.top().matrix,point_3);
                    print_triangle(scale(p1),scale(p2),scale(p3),1);



                }
                if(lineOfFile.compare("push")==0)
                {
                    push_pop_stack.push(matrixstack.top());
                    stack_size++;

                }
                if(lineOfFile.compare("pop")==0)
                {
                    if(stack_size==0) cout<<"Pop on Empty Stack"<<endl;
                    else
                    {

                        matrixstack.push(push_pop_stack.top());
                        push_pop_stack.pop();
                        stack_size--;
                    }
                }
                if(lineOfFile.compare("scale")==0)
                {

                    double point_x,point_y,point_z;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss(lineOfFile);
                    iss>>point_x>>point_y>>point_z;
                    struct point scale;
                    scale.x=point_x;
                    scale.y=point_y;
                    scale.z=point_z;
                    transMatrix tm_1;
                    tm_1.generate_scaling_matrix(scale);
                    transMatrix tm_2;
                    tm_2.matrix_multiplication(matrixstack.top().matrix,tm_1.matrix);
                    matrixstack.push(tm_2);

                }
                if(lineOfFile.compare("translate")==0)
                {

                    double point_x,point_y,point_z;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss(lineOfFile);
                    iss>>point_x>>point_y>>point_z;
                    struct point translate;
                    translate.x=point_x;
                    translate.y=point_y;
                    translate.z=point_z;
                    transMatrix tm_1;
                    tm_1.generate_translation_matrix(translate);
                    transMatrix tm_2;
                    tm_2.matrix_multiplication(matrixstack.top().matrix,tm_1.matrix);
                    matrixstack.push(tm_2);

                }
                if(lineOfFile.compare("rotate")==0)
                {

                    double angle,point_x,point_y,point_z;
                    getline(filein,lineOfFile);
                    linecount++;
                    std::istringstream iss(lineOfFile);
                    iss>>angle>>point_x>>point_y>>point_z;
                    struct point rota;
                    rota.x=point_x;
                    rota.y=point_y;
                    rota.z=point_z;
                    transMatrix tm_1;
                    tm_1.generate_rotation_matrix(rota,angle);
                    transMatrix tm_2;
                    tm_2.matrix_multiplication(matrixstack.top().matrix,tm_1.matrix);
                    matrixstack.push(tm_2);

                }
                if(lineOfFile.compare("end")==0)
                {

                    break;
                }

            }

        }
        filein.close();
        st_1.close();
        ifstream filein_2("stage_1.txt");
        st_2.open("stage_2.txt");
        if (!filein_2)
        {
            cout<<"Error in opening stage_1.txt"<<endl;
            exit(1);
        }
        else
        {
            point l,r,u;
            l=sub(look,eye);
            l=normalize(l);
            r=cross(l,up);
            r=normalize(r);
            u=cross(r,l);
            point negEye;
            negEye.x=-eyeX;
            negEye.y=-eyeY;
            negEye.z=-eyeZ;
            transMatrix tm_1,tm_2,tm_3;
            tm_1.generate_translation_matrix(negEye);
            tm_2.generate_rotation_matrix_for_view(l,r,u);
            tm_3.matrix_multiplication(tm_1.matrix,tm_2.matrix);
            int linecount=0;
            while (!filein_2.eof())
            {

                double point_x,point_y,point_z;
                getline(filein_2,lineOfFile);
                if(lineOfFile.empty()) break;

                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>point_x>>point_y>>point_z;
                struct point p1;
                p1.x=point_x;
                p1.y=point_y;
                p1.z=point_z;
                p1.w=1.0;
                getline(filein_2,lineOfFile);
                linecount++;
                std::istringstream iss_2(lineOfFile);
                iss_2>>point_x>>point_y>>point_z;
                struct point p2;
                p2.x=point_x;
                p2.y=point_y;
                p2.z=point_z;
                p2.w=1.0;
                getline(filein_2,lineOfFile);
                linecount++;
                std::istringstream iss_3(lineOfFile);
                iss_3>>point_x>>point_y>>point_z;
                struct point p3;
                p3.x=point_x;
                p3.y=point_y;
                p3.z=point_z;
                p3.w=1.0;
                point res1,res2,res3;
                res1=matrix_multi(tm_3.matrix,p1);
                res2=matrix_multi(tm_3.matrix,p2);
                res3=matrix_multi(tm_3.matrix,p3);
                print_triangle(scale(res1),scale(res2),scale(res3),2);
                getline(filein_2,lineOfFile);
                linecount++;

            }



        }
        filein_2.close();
        st_2.close();

        ifstream filein_3("stage_2.txt");
        st_3.open("stage_3.txt");
        if (!filein_3)
        {
            cout<<"Error in opening stage_2.txt"<<endl;
            exit(1);
        }
        else
        {
            double fov_x=fovY*aspectRatio;
            transMatrix tm1;
            //cout<<fov_x<<endl;
            tm1.generate_projection_matrix(fov_x);

            int linecount=0;
            while (!filein_3.eof())
            {

                double point_x,point_y,point_z;
                getline(filein_3,lineOfFile);
                if(lineOfFile.empty()) break;

                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>point_x>>point_y>>point_z;
                struct point p1;
                p1.x=point_x;
                p1.y=point_y;
                p1.z=point_z;
                p1.w=1.0;
                getline(filein_3,lineOfFile);
                linecount++;
                std::istringstream iss_2(lineOfFile);
                iss_2>>point_x>>point_y>>point_z;
                struct point p2;
                p2.x=point_x;
                p2.y=point_y;
                p2.z=point_z;
                p2.w=1.0;
                getline(filein_3,lineOfFile);
                linecount++;
                std::istringstream iss_3(lineOfFile);
                iss_3>>point_x>>point_y>>point_z;
                struct point p3;
                p3.x=point_x;
                p3.y=point_y;
                p3.z=point_z;
                p3.w=1.0;
                point res1,res2,res3;
                res1=matrix_multi(tm1.matrix,p1);
                res2=matrix_multi(tm1.matrix,p2);
                res3=matrix_multi(tm1.matrix,p3);
                print_triangle(scale(res1),scale(res2),scale(res3),3);
                getline(filein_3,lineOfFile);
                linecount++;

            }


        }
        filein_3.close();
        st_3.close();
        ifstream filein_4("stage_3.txt");
        //ifstream filein_4("temp.txt");
        if (!filein_4)
        {
            cout<<"Error in opening stage_3.txt"<<endl;
            exit(1);
        }
        else
        {
            int linecount=0;
            while (!filein_4.eof())
            {

                double point_x,point_y,point_z;
                getline(filein_4,lineOfFile);
                if(lineOfFile.empty()) break;

                linecount++;
                std::istringstream iss(lineOfFile);
                iss>>point_x>>point_y>>point_z;
                struct point p1;
                p1.x=point_x;
                p1.y=point_y;
                p1.z=point_z;
                p1.w=1.0;
                getline(filein_4,lineOfFile);
                linecount++;
                std::istringstream iss_2(lineOfFile);
                iss_2>>point_x>>point_y>>point_z;
                struct point p2;
                p2.x=point_x;
                p2.y=point_y;
                p2.z=point_z;
                p2.w=1.0;
                getline(filein_4,lineOfFile);
                linecount++;
                std::istringstream iss_3(lineOfFile);
                iss_3>>point_x>>point_y>>point_z;
                struct point p3;
                p3.x=point_x;
                p3.y=point_y;
                p3.z=point_z;
                p3.w=1.0;
                Triangle t;
                t.points[0]=p1;
                t.points[1]=p2;
                t.points[2]=p3;
                t.colors.r=rand()%256;
                t.colors.g=rand()%256;
                t.colors.b=rand()%256;
                triangles.push_back(t);
                getline(filein_4,lineOfFile);
                linecount++;

            }

        }
        filein_4.close();

        ifstream filein_5("config.txt");
        //st_3.open("stage_3.txt");
        if (!filein_5)
        {
            cout<<"Error in opening config.txt"<<endl;
            exit(1);
        }
        else
        {
            getline(filein_5,lineOfFile);
            std::istringstream iss(lineOfFile);
            iss>>screenWidth>>screenHeight;
            getline(filein_5,lineOfFile);
            std::istringstream iss2(lineOfFile);
            iss2>>left_lim_X;
            right_lim_X=-left_lim_X;

            getline(filein_5,lineOfFile);
            std::istringstream iss3(lineOfFile);
            iss3>>bottom_lim_Y;
            top_lim_Y=-bottom_lim_Y;
            getline(filein_5,lineOfFile);

            std::istringstream iss4(lineOfFile);
            iss4>>frontZ>>rearZ;

        }
        filein_5.close();

        double dx=(right_lim_X-left_lim_X)/screenWidth;
        double dy=(top_lim_Y-bottom_lim_Y)/screenHeight;
        double topY=top_lim_Y-(dy/2.0);
        double bottomY =-topY;
        double leftX=left_lim_X+(dx/2.0);
        double rightX=-leftX;

        double** ZBuffer=new double*[screenHeight];
        color** framebuffer=new color*[screenHeight];
        for(int i=0; i<screenHeight; i++)
        {
            ZBuffer[i]=new double[screenWidth];
        }
        for(int i=0; i<screenHeight; i++)
        {
            framebuffer[i]=new color[screenWidth];
        }
        for (int i=0; i<screenHeight; i++)
        {
            for(int j=0; j<screenWidth; j++)
            {
                ZBuffer[i][j]=rearZ;
            }
        }
        for (int i=0; i<screenHeight; i++)
        {
            for(int j=0; j<screenWidth; j++)
            {
                framebuffer[i][j].r=0;
                framebuffer[i][j].g=0;
                framebuffer[i][j].b=0;
            }
        }
        bitmap_image image(screenWidth,screenHeight);
        for (int i=0; i<screenWidth; i++)
        {
            for(int j=0; j<screenHeight; j++)
            {
                image.set_pixel(i,j,0,0,0);
            }
        }

        //cout<<triangle_vector.size()<<endl;

        for(int i=0; i<triangles.size(); i++)
        {
            int top_scan_line;
            int bottom_scan_line;
            double max_y,min_y;
            max_y = getMax(triangles[i].points[0].y,triangles[i].points[1].y, triangles[i].points[2].y);
            //cout<<minY<<"it is min"<<endl;
            int max_flag=0;
            if(max_y >= topY)
            {
                max_flag=1;
            }
            if(max_flag=1) top_scan_line=0;
            else top_scan_line=round((topY - max_y)/dy);
            min_y = getMin(triangles[i].points[0].y,triangles[i].points[1].y, triangles[i].points[2].y);
            int min_flag=0;
            if(min_y <= bottomY)
            {
                min_flag=1;

            }
            if(min_flag==1) bottom_scan_line = screenHeight - 1;
            else
            {
                int val=round((min_y +topY)/dy);
                bottom_scan_line = screenHeight-1-val;
            }

            for(int row=top_scan_line; row<=bottom_scan_line; row++)
            {

                double ys = topY - row*dy;
                point intersects[3];
                intersects[0].x=INF;
                intersects[0].y=ys;
                intersects[1].x=INF;
                intersects[1].y=ys;
                intersects[2].x=INF;
                intersects[2].y=ys;
                intersects[0].x=getintersection(triangles[i].points[0].x,triangles[i].points[0].y,triangles[i].points[1].x,triangles[i].points[1].y,ys);
                intersects[0].z=getintersection(triangles[i].points[0].z,triangles[i].points[0].y,triangles[i].points[1].z,triangles[i].points[1].y,ys);
                if(isInvalid(intersects[0],triangles[i].points[0].x,triangles[i].points[1].x,triangles[i].points[0].y,triangles[i].points[1].y))
                    intersects[0].x=INF;

                intersects[1].x=getintersection(triangles[i].points[1].x,triangles[i].points[1].y,triangles[i].points[2].x,triangles[i].points[2].y,ys);
                intersects[1].z=getintersection(triangles[i].points[1].z,triangles[i].points[1].y,triangles[i].points[2].z,triangles[i].points[2].y,ys);
                if(isInvalid(intersects[1],triangles[i].points[1].x,triangles[i].points[2].x,triangles[i].points[1].y,triangles[i].points[2].y))
                    intersects[1].x=INF;

                intersects[2].x=getintersection(triangles[i].points[2].x,triangles[i].points[2].y,triangles[i].points[0].x,triangles[i].points[0].y,ys);
                intersects[2].z=getintersection(triangles[i].points[2].z,triangles[i].points[2].y,triangles[i].points[0].z,triangles[i].points[0].y,ys);
                if(isInvalid(intersects[2],triangles[i].points[2].x,triangles[i].points[0].x,triangles[i].points[2].y,triangles[i].points[0].y))
                    intersects[2].x=INF;


                int max_ind, min_ind;
                max_ind=-1;
                min_ind=-1;

                double maxX=INF;
                double minX=-INF;
                int loop_flag=0;
                for(int s=0; s<3; s++)
                {
                    if(intersects[s].x!=INF)
                    {
                        max_ind=s;
                        min_ind=s;
                        maxX=intersects[s].x;
                        minX=intersects[s].x;
                        break;
                    }
                }
                for(int k=0; k<3; k++)
                {

                    if(intersects[k].x!=INF)
                    {
                        if(intersects[k].x<minX)
                        {
                            minX = intersects[k].x;
                            min_ind=k;

                        }
                        if(intersects[k].x>maxX)
                        {
                            max_ind=k;
                            maxX=intersects[k].x;

                        }
                    }
                }

                double za,zb;
                if(min_ind==0) za=intersects[0].z;
                else if(min_ind==1) za=intersects[1].z;
                else za=intersects[2].z;
                if(max_ind==0) zb=intersects[0].z;
                else if(max_ind==1) zb=intersects[1].z;
                else zb=intersects[2].z;
                int leftcol,rightcol;
                double maxVal=intersects[max_ind].x;
                double minVal=intersects[min_ind].x;
                int min_flag=0;
                if(minVal<=leftX)
                {
                    min_flag=1;
                }
                if (min_flag==1) leftcol=0;
                else leftcol=round((minVal+rightX)/dx);
                int max_flag=0;
                if(maxVal>=rightX)
                {
                    max_flag=1;

                }
                if(max_flag==1) rightcol=screenWidth-1;
                else
                {
                    int val=round((rightX-maxVal)/dx);
                    rightcol=screenWidth-1-val;
                }
                double xa=minX;
                double xb=maxX;


                double cons = dx*(zb - za)/(xb - xa);
                double zp;
                int colflag=0;
                for(int col=leftcol; col<=rightcol; col++)
                {
                    double xp=leftX + leftcol*dx;
                    if(colflag==0)
                    {
                        zp=zb-(xb-xp)*(zb - za)/(xb- xa);
                        colflag=1;
                    }
                    else
                    {
                        zp=zp+cons;
                    }
                    if(zp>frontZ &&zp<ZBuffer[row][col])
                    {
                        framebuffer[row][col].r=triangles[i].colors.r;
                        framebuffer[row][col].g=triangles[i].colors.g;
                        framebuffer[row][col].b=triangles[i].colors.b;
                        ZBuffer[row][col]=zp;

                    }
                }
            }
        }
        for(int i=0; i<screenHeight; i++)
        {
            for(int j=0; j<screenWidth; j++)
            {
                image.set_pixel(j,i,framebuffer[i][j].r,framebuffer[i][j].g,framebuffer[i][j].b);
            }
        }
        image.save_image("out.bmp");

        st_4.open("z_buffer.txt");
        for(int i=0; i<screenHeight; i++)
        {
            for(int j=0; j<screenWidth; j++)
            {
                if(ZBuffer[i][j]<rearZ)
                    st_4<<fixed<<setprecision(6)<<ZBuffer[i][j]<<"\t";

            }
            st_4<<endl;
        }
        st_4.close();

        for(int i=0; i<screenHeight; i++) delete[] ZBuffer[i];
        delete[] ZBuffer;
        for(int i=0; i<screenHeight; i++) delete[] framebuffer[i];
        delete[] framebuffer;
        vector<Triangle>().swap(triangles);

    }

    return 0;
}



