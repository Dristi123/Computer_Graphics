#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include<bits/stdc++.h>
#include<iostream>
#define INF numeric_limits<double>::infinity()
#define pi (2*acos(0.0))
#include<stdio.h>
#include<stdlib.h>

using namespace std;

struct Vector3D
{
    double x;
    double y;
    double z;

    Vector3D(double a,double b,double c)
    {
        this->x=a;
        this->y=b;
        this->z=c;
    }
    Vector3D()
    {

    }
};
struct Vector3D normalize(struct Vector3D a)
{
    Vector3D temp;
    double r=sqrt((a.x*a.x)+(a.y*a.y)+(a.z*a.z));
    temp.x=a.x/r;
    temp.y=a.y/r;
    temp.z=a.z/r;
    return temp;
}

double dot(Vector3D a,Vector3D b)
{
    double p=(a.x*b.x)+(a.y*b.y)+(a.z*b.z);
    return p;
}

void initializecolor(double* colarr) {
   for(int i=0;i<3;i++) colarr[i]=0.0;
   }
class Light
{
public:
    Vector3D light_position;
    double color[3];
    Vector3D direction;
    double cut_off_angle;
    Light() {}
    Light(double pos_x,double pos_y,double pos_z)
    {
        Vector3D temp;
        temp.x=pos_x;
        temp.y=pos_y;
        temp.z=pos_z;
        this->light_position=temp;
    }
    void setColor(double r,double g,double b)
    {
        this->color[0]=r;
        this->color[1]=g;
        this->color[2]=b;
    }
    void setDirection(Vector3D dir)
    {
        this->direction=dir;
    }
    void setAngle(double angle)
    {
        this->cut_off_angle=angle;
    }
    void draw()
    {
        //cout<<cut_off_angle<<endl;
        if(cut_off_angle==INF)
        {
            glPointSize(7);
            glColor3f(color[0],color[1],color[2]);
            glBegin(GL_POINTS);
            {
                glVertex3f(light_position.x,light_position.y,light_position.z);
            }
            glEnd();
        }
        else
        {
            //glPointSize(10);
            glColor3f(color[0],color[1],color[2]);
            //glBegin(GL_POINTS);
            //{
            //glVertex3f(light_position.x,light_position.y,light_position.z);
            //}
            glBegin(GL_QUADS);
            {
                glVertex3f( light_position.x, light_position.y,light_position.z);
                glVertex3f( light_position.x+15,light_position.y,light_position.z);
                glVertex3f(light_position.x+15,light_position.y+15,light_position.z);
                glVertex3f(light_position.x, light_position.y+15,light_position.z);
            }
            glEnd();
        }

    }
};
class Ray
{
public:
    Vector3D Ro;
    Vector3D Rd;

    Ray(Vector3D a,Vector3D b)
    {

        Setter(a,b);
    }
    Ray()
    {

    }
    void Setter(Vector3D a, Vector3D b)
    {
    this->Ro=a;
    this->Rd=b;
    }
};
Vector3D multiplyScalar(Vector3D a, double b)
{
    Vector3D temp;
    temp.x=a.x*b;
    temp.y=a.y*b;
    temp.z=a.z*b;
    return temp;
}
Vector3D createNewPoint_Add(Vector3D a, Vector3D b)
{
    Vector3D temp;
    temp.x=a.x+b.x;
    temp.y=a.y+b.y;
    temp.z=a.z+b.z;
    return temp;
}
Vector3D createNewPoint_Sub(Vector3D a, Vector3D b)
{
    Vector3D temp;
    temp.x=a.x-b.x;
    temp.y=a.y-b.y;
    temp.z=a.z-b.z;
    return temp;
}
void printPoint(Vector3D point)
{
    cout<<point.x<<endl;
    cout<<point.y<<endl;
    cout<<point.z<<endl;
}
Vector3D checkmax(Vector3D a,Vector3D b) {
    if(a.x>b.x)return a;
    else return b;
    }
Vector3D reflected_ray(Ray *r,struct Vector3D normal)
{
    Vector3D refl;
    double dot_prod=(dot(r->Rd,normal)*2.0);
    Vector3D inter=multiplyScalar(normal,dot_prod);
    refl=createNewPoint_Sub(r->Rd,inter);
    return normalize(refl);
}
void delete_array(double *arr) {
 delete[] arr;
 }
double check_angle(Vector3D a,Vector3D b)
{
    double mag_1=sqrt(pow(a.x,2)+pow(a.y,2)+pow(a.z,2));
    double mag_2=sqrt(pow(b.x,2)+pow(b.y,2)+pow(b.z,2));
    a=normalize(a);
    b=normalize(b);
    double angle=acos(dot(a,b));
    angle=angle*(180.0/pi);
    return angle;
}
Vector3D get_light_pos(Vector3D dir,Vector3D inter_point)
{
    return createNewPoint_Add(inter_point,multiplyScalar(dir,0.000001));
}
class Object
{
public:
    Vector3D reference_point;
    double height;
    double width;
    double length;
    int shine;
    double color[3];
    double coEfficients[4];
    Object()
    {

    }
    virtual void draw()
    {

    }
    void setCoEfficients(double x,double y,double z,double w)
    {
        this->coEfficients[0]=x;
        this->coEfficients[1]=y;
        this->coEfficients[2]=z;
        this->coEfficients[3]=w;
    }
    void setColor(double red,double green,double blue)
    {
        this->color[0]=red;
        this->color[1]=green;
        this->color[2]=blue;
    }
    void setShine(int sh)
    {
        this->shine=sh;
    }
    virtual double intersect(Ray *r, double *color_pointer_sphere, int level)
    {
        return -1.0;
    }
    virtual double get_parameter_t(Ray *r)
    {

    }
    virtual void computePhong_and_lambert(Ray inc,Vector3D normal,Ray *r, Light *light,double* color_pointer_sphere) {

    }

};

int level_of_rec;
vector<Light*> lights;
vector<Object*> objects;
struct Vector3D pos;
Vector3D cross(Vector3D a,Vector3D b)
{
    Vector3D result;
    result.x=(a.y*b.z)-(b.y*a.z);
    result.y=(a.z*b.x)-(a.x*b.z);
    result.z=(a.x*b.y)-(a.y*b.x);
    return result;
}
bool check_t(double t,double min_t)
{
    if(t>0 && t<min_t) return true;
    else return false;
}
class Sphere:public Object
{
private:
    int slices;
    int stacks;
public:
    Sphere(Vector3D centre,double radius)
    {
        reference_point=centre;
        length=radius;
        this->slices=50;
        this->stacks=50;
    }
    void computeAmbianceColor(double* color_pointer_sphere)
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_sphere[i]=this->color[i]*this->coEfficients[0];
        }
    }
    bool check_sphere_normal(Vector3D normal) {
        if(normal.z==0) return true;
    }
    void computeReflectionColor(double* color_pointer_sphere, Light* light, double lambet,double phong)
    {
        double temp_1,temp_2,temp_3;
        double lambert_component=max(lambet,0.0)*this->coEfficients[1];
        double phong_component=pow(max(phong,0.0),this->shine)*this->coEfficients[2];
        temp_1=(light->color[0]*lambert_component*this->color[0])+(light->color[0]*phong_component*this->color[0]);
        temp_2=(light->color[1]*lambert_component*this->color[1])+(light->color[1]*phong_component*this->color[1]);
        temp_3=(light->color[2]*lambert_component*this->color[2])+(light->color[2]*phong_component*this->color[2]);
        color_pointer_sphere[0]=color_pointer_sphere[0]+temp_1;
        color_pointer_sphere[1]=color_pointer_sphere[1]+temp_2;
        color_pointer_sphere[2]=color_pointer_sphere[2]+temp_3;
    }
    void computeRecursiveColor(double* color_pointer_sphere,double* ref_col )
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_sphere[i]=color_pointer_sphere[i]+ref_col[i]*this->coEfficients[3];
        }
    }

    double get_parameter_t(Ray *r)
    {
        double a,b,c,t,dis;
        Vector3D l;
        l.x=r->Ro.x-reference_point.x;
        l.y=r->Ro.y-reference_point.y;
        l.z=r->Ro.z-reference_point.z;
        c=(dot(l,l)-length*length);
        a=dot(r->Rd,r->Rd);
        b=(dot(r->Ro,r->Rd)-dot(r->Rd,reference_point));
        b=b*2.0;
        dis=(b*b)-(4*a*c*1.0);
        if(dis<0)
        {
            t=-1;
        }
        else if(dis>0)
        {
            double root1,root2;
            double temp=-b;
            double tem_nume=(2.0*a);
            root1=temp+sqrt(dis);
            root1=root1/tem_nume;
            //root1=(-b+sqrt(dis))/(2.0*a);
            root2=temp-sqrt(dis);
            root2=root2/tem_nume;
            if(root1>0)
            {
                if(root2>0)
                {
                    t= min(root1,root2);
                }
                else t=root1;
            }
            else
            {
                if(root2>0) t=root2;
                else t=-1;
            }

        }
        else
        {
            t=-b/(2.0*a);
        }
        return t;

    }
    void computePhong_and_lambert(Ray inc,Vector3D normal,Ray *r, Light *light,double* color_pointer_sphere)
    {
        double lambert=dot(inc.Rd,normal);
        Vector3D reflected=reflected_ray(&inc,normal);
        double phong=dot(reflected,r->Rd);
        this->computeReflectionColor(color_pointer_sphere,light,lambert,phong);
    }
    void draw()
    {

        Vector3D points[100][100];
        int i,j;
        double h,r;

        //generate points
        for(i=0; i<=stacks; i++)
        {
            h=this->length*sin(((double)i/(double)stacks)*(pi/2));
            r=this->length*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0; j<=slices; j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }
        //draw quads using generated points
        for(i=0; i<stacks; i++)
        {
            //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
            //glColor3f(1,0,0);
            glColor3f(color[0], color[1], color[2]);
            for(j=0; j<slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    //upper hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,points[i][j].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,points[i+1][j].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,points[i][j+1].z+reference_point.z);

                    //lower hemisphere
                    glVertex3f(points[i][j].x+reference_point.x,points[i][j].y+reference_point.y,-points[i][j].z+reference_point.z);
                    glVertex3f(points[i+1][j].x+reference_point.x,points[i+1][j].y+reference_point.y,-points[i+1][j].z+reference_point.z);
                    glVertex3f(points[i+1][j+1].x+reference_point.x,points[i+1][j+1].y+reference_point.y,-points[i+1][j+1].z+reference_point.z);
                    glVertex3f(points[i][j+1].x+reference_point.x,points[i][j+1].y+reference_point.y,-points[i][j+1].z+reference_point.z);

                }
                glEnd();
            }
        }
    }


    double intersect(Ray *,double *,int);

};
int check_if_obscured(Ray r,int length) {
    int flag=0;
     for(int i=0;i<objects.size();i++) {
        double param_t=objects[i]->get_parameter_t(&r);
                if(param_t<=0) continue;
                else
                {
                    if(param_t<length)
                    {
                        flag=1;
                        break;
                    }
                }
     }
     return flag;
}
double Sphere::intersect(Ray *r, double *color_pointer_sphere, int level)
{

    double t;
    t=get_parameter_t(r);
    if(t>0)
    {
        if(level==0)
        {
            return t;
        }

        Vector3D intersection_point((r->Ro.x+(r->Rd.x*t)),(r->Ro.y+(r->Rd.y*t)),(r->Ro.z+(r->Rd.z*t)));
        double norm_x,norm_y,norm_z;
        Vector3D normal=createNewPoint_Sub(intersection_point,reference_point);
        normal=normalize(normal);
        double distance=sqrt(pow(pos.x-reference_point.x, 2.0)+pow(pos.y-reference_point.y, 2.0)+pow(pos.z-reference_point.z, 2.0));
        if(distance<=length)
        {
            normal.x=normal.x*(-1.0);
            normal.y=normal.y*(-1.0);
            normal.z=normal.z*(-1.0);
        }
        computeAmbianceColor(color_pointer_sphere);
        for(int i=0; i<lights.size(); i++)
        {

            Vector3D light_dir=createNewPoint_Sub(lights[i]->light_position,intersection_point);
            Vector3D temp=light_dir;

            light_dir=normalize(light_dir);
            if(lights[i]->cut_off_angle!=INF)
            {
                Vector3D temp=multiplyScalar(light_dir,-1.0);
                double angle=check_angle(temp,lights[i]->direction);
                if(angle>lights[i]->cut_off_angle) continue;
            }
            Vector3D light_pos=get_light_pos(light_dir,intersection_point);
            Ray incident_ray(light_pos,light_dir);

            double length_of_light_dir=sqrt((temp.x*temp.x)+(temp.y*temp.y)+(temp.z*temp.z));
            int flag=check_if_obscured(incident_ray,length_of_light_dir);
            //check if ray is obscured
            if(flag==0)
            {
                 double lambert=dot(incident_ray.Rd,normal);
                 Vector3D reflected=reflected_ray(&incident_ray,normal);
                 double phong=dot(reflected,r->Rd);
                 computeReflectionColor(color_pointer_sphere,lights[i],lambert,phong);

            }

        }

        if(level>=level_of_rec) return t;

        else
        {

            Vector3D refl_again=reflected_ray(r,normal);
            Vector3D refl_pos=get_light_pos(refl_again,intersection_point);
            Ray Ref_ray(refl_pos,refl_again);
            double tMin=INF;
            double t;
            int nearest=-9999;
            for(int k=0; k<objects.size(); k++)
            {
                double *ref_Color=new double[3];
                initializecolor(ref_Color);
                t=objects[k]->intersect(&Ref_ray,ref_Color,0);
                if(check_t(t,tMin))
                {

                    nearest=k;
                    tMin=t;

                }
                delete_array(ref_Color);
            }
            if(nearest!=-9999)
            {
                double *ref_Color=new double[3];
                initializecolor(ref_Color);
                tMin=objects[nearest]->intersect(&Ref_ray,ref_Color,level+1);
                if(tMin!=-1)
                {
                    computeRecursiveColor(color_pointer_sphere,ref_Color);
                }
                delete_array(ref_Color);

            }
            color_pointer_sphere[0]>1?color_pointer_sphere[0]=1:(color_pointer_sphere[0]<0?color_pointer_sphere[0]=0:color_pointer_sphere[0]=color_pointer_sphere[0]);
            color_pointer_sphere[1]>1?color_pointer_sphere[1]=1:(color_pointer_sphere[1]<0?color_pointer_sphere[1]=0:color_pointer_sphere[1]=color_pointer_sphere[1]);
            color_pointer_sphere[2]>1?color_pointer_sphere[2]=1:(color_pointer_sphere[2]<0?color_pointer_sphere[2]=0:color_pointer_sphere[2]=color_pointer_sphere[2]);
            return tMin;
        }
    }
    else return -1;

}
class Triangle:public Object
{
private:
    Vector3D side1;
    Vector3D side2;
    Vector3D side3;
public:
    Vector3D vertex_1;
    Vector3D vertex_2;
    Vector3D vertex_3;

    Triangle(double a_x,double a_y,double a_z,double b_x,double b_y,double b_z,double c_x,double c_y,double c_z)
    {
        vertex_1=Vector3D(a_x,a_y,a_z);
        vertex_2=Vector3D(b_x,b_y,b_z);
        vertex_3=Vector3D(c_x,c_y,c_z);
        this->reference_point=Vector3D(0.0,0.0,0.0);
        this->length=0.0;
        this->side1=createNewPoint_Sub(vertex_1,vertex_2);
        this->side2=createNewPoint_Sub(vertex_1,vertex_3);
        this->side3=createNewPoint_Sub(vertex_2,vertex_3);
    }
    void draw()
    {

        glBegin(GL_TRIANGLES);
        {
            glColor3f(color[0],color[1],color[2]);
            glVertex3f(vertex_1.x,vertex_1.y,vertex_1.z);
            glVertex3f(vertex_2.x,vertex_2.y,vertex_2.z);
            glVertex3f(vertex_3.x,vertex_3.y,vertex_3.z);
        }
        glEnd();
    }
    Vector3D get_cross_of_two_sides()
    {
        Vector3D edge_1=multiplyScalar(this->side1,-1.0);
        Vector3D edge_2=multiplyScalar(this->side2,-1.0);
        return normalize(cross(edge_1,edge_2));
    }
    bool check_normal(Vector3D intersecting_point,Vector3D normal)
    {
        if(dot(intersecting_point,normal)>0) return true;
        else return false;
    }
    void computeAmbianceColor(double* color_pointer)
    {
        for(int i=0; i<3; i++)
        {
            color_pointer[i]=this->color[i]*this->coEfficients[0];
        }
    }
    void computePhong_and_lambert(Ray inc,Vector3D normal,Ray *r, Light *light,double* color_pointer_sphere)
    {
        double lambert=dot(inc.Rd,normal);
        Vector3D reflected=reflected_ray(&inc,normal);
        double phong=dot(reflected,r->Rd);
        this->computeReflectionColor(color_pointer_sphere,light,lambert,phong);
    }
    void computeReflectionColor(double* color_pointer, Light* light, double lambet,double phong)
    {

        double temp_1,temp_2,temp_3;
        double lambert_component=max(lambet,0.0)*this->coEfficients[1];
        double phong_component=pow(max(phong,0.0),this->shine)*this->coEfficients[2];
        temp_1=(light->color[0]*lambert_component*this->color[0])+(light->color[0]*phong_component*this->color[0]);
        temp_2=(light->color[1]*lambert_component*this->color[1])+(light->color[1]*phong_component*this->color[1]);
        temp_3=(light->color[2]*lambert_component*this->color[2])+(light->color[2]*phong_component*this->color[2]);
        color_pointer[0]=color_pointer[0]+temp_1;
        color_pointer[1]=color_pointer[1]+temp_2;
        color_pointer[2]=color_pointer[2]+temp_3;
    }
    void computeRecursiveColor(double* color_pointer,double* ref_col )
    {
        for(int i=0; i<3; i++)
        {
            color_pointer[i]=color_pointer[i]+ref_col[i]*this->coEfficients[3];
        }
    }
    bool check_triangle_norm(Vector3D norm) {
    if(norm.z==0) return true;
    }
    double get_parameter_t(Ray *r)
    {

        Vector3D co_eff_1,co_eff_2,co_eff_3;
        co_eff_1=this->side1;
        co_eff_2=this->side2;
        co_eff_3=createNewPoint_Sub(vertex_1,r->Ro);
        double determininant;
        double rdx=r->Rd.x;
        double rdy=r->Rd.y;
        double rdz=r->Rd.z;
        determininant=(co_eff_1.x*((co_eff_2.y*rdz)-(co_eff_2.z*rdy)))+(co_eff_1.y*((co_eff_2.z*rdx)-(co_eff_2.x*rdz)))+(co_eff_1.z*((co_eff_2.x*rdy)-(co_eff_2.y*rdx)));

        double beta,T,gamma;
        beta=co_eff_3.x*((co_eff_2.y*r->Rd.z)-(co_eff_2.z*r->Rd.y));
        beta=beta+(co_eff_2.x*((co_eff_3.z*r->Rd.y)-(co_eff_3.y*r->Rd.z)));
        beta=beta+(r->Rd.x*((co_eff_3.y*co_eff_2.z)-(co_eff_3.z*co_eff_2.y)));
        beta=beta/determininant;


        gamma=co_eff_1.x*((co_eff_3.y*r->Rd.z)-(co_eff_3.z*r->Rd.y));
        gamma=gamma+(co_eff_3.x*((co_eff_1.z*r->Rd.y)-(co_eff_1.y*r->Rd.z)));
        gamma=gamma+(r->Rd.x*((co_eff_1.y*co_eff_3.z)-(co_eff_1.z*co_eff_3.y)));
        gamma=gamma/determininant;
        //cout<<Beta<<endl;
        //cout<<"Beta"<<Beta<<endl;
        //cout<<"amar"<<beta<<endl;
        //cout<<"Gamma"<<Gamma<<endl;
        //cout<<"amar"<<gamma<<endl;
        T=co_eff_1.x*((co_eff_2.y*co_eff_3.z)-(co_eff_3.y*co_eff_2.z));
        T=T+(co_eff_2.x*((co_eff_1.z*co_eff_3.y)-(co_eff_1.y*co_eff_3.z)));
        T=T+(co_eff_3.x*((co_eff_1.y*co_eff_2.z)-(co_eff_1.z*co_eff_2.y)));
        T/=determininant;
        if(determininant==0)
        {
            return -1;
        }
        else
        {
            if(beta>0.0 && gamma>0.0 && (beta+gamma<1.0))
            {
                if(T<0) return -1;
                else return T;
            }
            else
            {

                return -1;
            }
        }
        //cout<<"T"<<t<<endl;
        //cout<<"amar"<<T<<endl;
        //cout<<t<<endl;



    }
    double intersect(Ray *,double *,int);
};
double Triangle::intersect(Ray *r, double *color_pointer_triangle, int level)
{
    double t;
    t=get_parameter_t(r);
    if(t>0)
    {

        //cout<<t<<endl;
        if(level==0)
        {
            return t;
        }

        this->computeAmbianceColor(color_pointer_triangle);
        Vector3D intersection_point((r->Ro.x+(r->Rd.x*t)),(r->Ro.y+(r->Rd.y*t)),(r->Ro.z+(r->Rd.z*t)));
        Vector3D normal;
        normal=this->get_cross_of_two_sides();
        Vector3D temp=multiplyScalar(r->Rd,-1.0);
        if(dot(temp,normal)<0.0) normal=multiplyScalar(normal,-1.0);
        //double distance=sqrt(pow(pos.x-reference_point.x, 2.0)+pow(pos.y-reference_point.y, 2.0)+pow(pos.z-reference_point.z, 2.0));
        //if (check_normal(intersection_point,normal))
        //{
            //normal=multiplyScalar(normal,-1.0);
        //}

        for(int i=0; i<lights.size(); i++)
        {

            Vector3D light_dir=createNewPoint_Sub(lights[i]->light_position,intersection_point);
            Vector3D temp=light_dir;

            light_dir=normalize(light_dir);
            if(lights[i]->cut_off_angle!=INF)
            {
                Vector3D temp=multiplyScalar(light_dir,-1.0);
                double angle=check_angle(temp,lights[i]->direction);
                if(angle>lights[i]->cut_off_angle) continue;
            }
            Vector3D light_pos=get_light_pos(light_dir,intersection_point);
            Ray incident_ray(light_pos,light_dir);
            //int flag=0;
            double length_of_light_dir=sqrt((temp.x*temp.x)+(temp.y*temp.y)+(temp.z*temp.z));
            //check if ray is obscured

            int flag=check_if_obscured(incident_ray,length_of_light_dir);
            if(flag==0)
            {
                double lambert=dot(incident_ray.Rd,normal);
                 Vector3D reflected=reflected_ray(&incident_ray,normal);
                 double phong=dot(reflected,r->Rd);
                 computeReflectionColor(color_pointer_triangle,lights[i],lambert,phong);
            }

        }

        if(level>=level_of_rec) return t;

        else
        {

            Vector3D refl_again=reflected_ray(r,normal);
            Vector3D refl_pos=get_light_pos(refl_again,intersection_point);
            Ray Ref_ray(refl_pos,refl_again);
            double tMin=INF;
            double t;
            int nearest=-9999;
            for(int k=0; k<objects.size(); k++)
            {
                double *ref_Color=new double[3];
                ref_Color[0]=0.0;
                ref_Color[1]=0.0;
                ref_Color[2]=0.0;
                t=objects[k]->intersect(&Ref_ray,ref_Color,0);
                if(check_t(t,tMin))
                {

                    nearest=k;
                    tMin=t;

                }
                delete_array(ref_Color);
            }
            if(nearest!=-9999)
            {
                double *ref_Color=new double[3];
                ref_Color[0]=0.0;
                ref_Color[1]=0.0;
                ref_Color[2]=0.0;
                tMin=objects[nearest]->intersect(&Ref_ray,ref_Color,level+1);
                if(tMin!=-1)
                {
                    computeRecursiveColor(color_pointer_triangle,ref_Color);
                }
                delete_array(ref_Color);

            }
            color_pointer_triangle[0]>1?color_pointer_triangle[0]=1:(color_pointer_triangle[0]<0?color_pointer_triangle[0]=0:color_pointer_triangle[0]=color_pointer_triangle[0]);
            color_pointer_triangle[1]>1?color_pointer_triangle[1]=1:(color_pointer_triangle[1]<0?color_pointer_triangle[1]=0:color_pointer_triangle[1]=color_pointer_triangle[1]);
            color_pointer_triangle[2]>1?color_pointer_triangle[2]=1:(color_pointer_triangle[2]<0?color_pointer_triangle[2]=0:color_pointer_triangle[2]=color_pointer_triangle[2]);
            return tMin;
        }
    }
    else return -1;
}
class Floor:public Object
{
private:
    double florr_width;
    double tile_width;
public:

    Floor(double floor_width,double tile_width)
    {
        length=tile_width;
        double val=floor_width/2.0;
        val=(-1.0)*val;
        Vector3D temp_point(val,val,0);
        reference_point=temp_point;
        this->florr_width=floor_width;
        this->tile_width=tile_width;
        setColor(0.0,0.0,0.0);

    }
    void computeAmbianceColor(double* color_pointer_sphere)
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_sphere[i]=this->color[i]*this->coEfficients[0];
        }
    }
    void computeReflectionColor(double* color_pointer_sphere, Light* light, double lambet,double phong)
    {
        double temp_1,temp_2,temp_3;
        double lambert_component=max(lambet,0.0)*this->coEfficients[1];
        double phong_component=pow(max(phong,0.0),this->shine)*this->coEfficients[2];
        temp_1=(light->color[0]*lambert_component*this->color[0])+(light->color[0]*phong_component*this->color[0]);
        temp_2=(light->color[1]*lambert_component*this->color[1])+(light->color[1]*phong_component*this->color[1]);
        temp_3=(light->color[2]*lambert_component*this->color[2])+(light->color[2]*phong_component*this->color[2]);
        color_pointer_sphere[0]=color_pointer_sphere[0]+temp_1;
        color_pointer_sphere[1]=color_pointer_sphere[1]+temp_2;
        color_pointer_sphere[2]=color_pointer_sphere[2]+temp_3;
    }
    void computePhong_and_lambert(Ray inc,Vector3D normal,Ray *r, Light *light,double* color_pointer_sphere)
    {
        double lambert=dot(inc.Rd,normal);
        Vector3D reflected=reflected_ray(&inc,normal);
        double phong=dot(reflected,r->Rd);
        this->computeReflectionColor(color_pointer_sphere,light,lambert,phong);
    }
    void computeRecursiveColor(double* color_pointer_sphere,double* ref_col )
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_sphere[i]=color_pointer_sphere[i]+ref_col[i]*this->coEfficients[3];
        }
    }
    bool if_normal_valid(Vector3D nor) {
    if(nor.z==0) return true;
    }
    Vector3D check_normal(Vector3D normal)
    {
        if(dot(pos,normal)>0) return normal;
        else return multiplyScalar(normal,-1.0);
    }
    void draw()
    {

        bool col_flag = false;
        int no_of_rows=(int)(florr_width/tile_width);
        double val=-florr_width*0.5;
        for(int i=0; i<no_of_rows; i++)
        {
            for(int j=0; j<no_of_rows; j++)
            {
                Vector3D temp(val+tile_width*j, val+tile_width*i, 0.0);
                if(col_flag==true)
                {
                    glColor3f(0, 0, 0);
                }
                else
                {
                    glColor3f(1, 1, 1);
                }


                glBegin(GL_QUADS);
                {
                    glVertex3f(temp.x, temp.y, temp.z);
                    glVertex3f(temp.x+length, temp.y, temp.z);
                    glVertex3f(temp.x + length, temp.y + length, temp.z);
                    glVertex3f(temp.x, temp.y+length, temp.z);
                }
                glEnd();
                col_flag=!col_flag;

            }

            col_flag=!col_flag;
        }
    }
    bool isIntersectionPointInavlid(Vector3D point)
    {
        if(point.x>-reference_point.x) return true;
        if(point.x<reference_point.x)  return true;
        if(point.y>-reference_point.y) return true;
        if(point.y<reference_point.y)  return true;
        return false;
    }

    double get_parameter_t(Ray* r)
    {
        Vector3D normal(0.0,0.0,1.0);
        double temp=dot(normal,r->Ro)/dot(normal,r->Rd);
        return (-1.0)*temp;

    }
    double intersect(Ray *,double *,int);

};
double Floor::intersect(Ray *r, double *color_pointer_floor, int level)
{
    Vector3D normal(0.0,0.0,1.0);
    normal=check_normal(normal);
    if(r->Rd.z==0) return -1;
    double t;
    t=get_parameter_t(r);
    if(t>0)
    {


        Vector3D intersection_point((r->Ro.x+(r->Rd.x*t)),(r->Ro.y+(r->Rd.y*t)),(r->Ro.z+(r->Rd.z*t)));
        if(isIntersectionPointInavlid(intersection_point))
        {
            return -1;
        }
        if(level==0)
        {
            return t;
        }
        Vector3D temp_ref_point=createNewPoint_Sub(intersection_point,reference_point);
        int a,b;
        a = (int)(temp_ref_point.x/length);
        b = (int)(temp_ref_point.y/length);
        double val=(((a+b)%2==0)?1.0:0.0);
        this->color[0]=val;
        this->color[1]=this->color[0];
        this->color[2]=this->color[1];

        computeAmbianceColor(color_pointer_floor);

        for(int i=0; i<lights.size(); i++)
        {

            Vector3D light_dir=createNewPoint_Sub(lights[i]->light_position,intersection_point);
            Vector3D temp=light_dir;
            light_dir=normalize(light_dir);
            if(lights[i]->cut_off_angle!=INF)
            {
                Vector3D temp=multiplyScalar(light_dir,-1.0);
                double angle=check_angle(temp,lights[i]->direction);
                if(angle>lights[i]->cut_off_angle) continue;
            }
            Vector3D light_pos=get_light_pos(light_dir,intersection_point);
            Ray incident_ray(light_pos,light_dir);

            double length_of_light_dir=sqrt((temp.x*temp.x)+(temp.y*temp.y)+(temp.z*temp.z));
            //check if ray is obscured
            int flag=check_if_obscured(incident_ray,length_of_light_dir);
            if(flag==0)
            {
                double lambert=dot(incident_ray.Rd,normal);
                 Vector3D reflected=reflected_ray(&incident_ray,normal);
                 double phong=dot(reflected,r->Rd);
                 computeReflectionColor(color_pointer_floor,lights[i],lambert,phong);
            }

        }

        if(level>=level_of_rec) return t;

        else
        {

            Vector3D refl_again=reflected_ray(r,normal);
            Vector3D refl_pos=get_light_pos(refl_again,intersection_point);
            Ray Ref_ray(refl_pos,refl_again);
            double tMin=INF;
            double t;
            int nearest=-9999;
            for(int k=0; k<objects.size(); k++)
            {
                double *ref_Color=new double[3];
                initializecolor(ref_Color);
                t=objects[k]->intersect(&Ref_ray,ref_Color,0);
                if(check_t(t,tMin))
                {

                    nearest=k;
                    tMin=t;

                }
                delete_array(ref_Color);
            }
            if(nearest!=-9999)
            {
                double *ref_Color=new double[3];
                initializecolor(ref_Color);
                tMin=objects[nearest]->intersect(&Ref_ray,ref_Color,level+1);
                if(tMin!=-1)
                {
                    computeRecursiveColor(color_pointer_floor,ref_Color);
                }
                delete_array(ref_Color);
            }
            color_pointer_floor[0]>1?color_pointer_floor[0]=1:(color_pointer_floor[0]<0?color_pointer_floor[0]=0:color_pointer_floor[0]=color_pointer_floor[0]);
            color_pointer_floor[1]>1?color_pointer_floor[1]=1:(color_pointer_floor[1]<0?color_pointer_floor[1]=0:color_pointer_floor[1]=color_pointer_floor[1]);
            color_pointer_floor[2]>1?color_pointer_floor[2]=1:(color_pointer_floor[2]<0?color_pointer_floor[2]=0:color_pointer_floor[2]=color_pointer_floor[2]);
            return tMin;
        }
    }
    else return -1;
}
class GenCoeffs
{
private:
    double a;
    double b;
    double c;
    double d;
    double e;
    double f;
    double g;
    double h;
    double i;
    double j;

public:
    GenCoeffs()
    {
    }
    GenCoeffs(double a,double b,double c,double d,double e,double f,double g,double h,double i,double j)
    {
        this->a=a;
        this->b=b;
        this->c=c;
        this->d=d;
        this->e=e;
        this->f=f;
        this->g=g;
        this->h=h;
        this->i=i;
        this->j=j;
    }
    double get_a()
    {
        return a;
    }
    double get_b()
    {
        return b;
    }
    double get_c()
    {
        return c;
    }
    double get_d()
    {
        return d;
    }
    double get_e()
    {
        return e;
    }
    double get_f()
    {
        return f;
    }
    double get_g()
    {
        return g;
    }
    double get_h()
    {
        return h;
    }
    double get_i()
    {
        return i;
    }
    double get_j()
    {
        return j;
    }
};
class General_quadratic:public Object
{
public:

    GenCoeffs parameters;
    General_quadratic(GenCoeffs parameters,double x,double y,double z,double length,double width,double height)
    {
        this->parameters=parameters;
        this->length=length;
        this->width=width;
        this->height=height;
        Vector3D reference_cube(x,y,z);
        this->reference_point=reference_cube;
    }
    void computeAmbianceColor(double* color_pointer_quad)
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_quad[i]=this->color[i]*this->coEfficients[0];
        }
    }
    void computeReflectionColor(double* color_pointer_quad, Light* light, double lambet,double phong)
    {
        double temp_1,temp_2,temp_3;
        double lambert_component=max(lambet,0.0)*this->coEfficients[1];
        double phong_component=pow(max(phong,0.0),this->shine)*this->coEfficients[2];
        temp_1=(light->color[0]*lambert_component*this->color[0])+(light->color[0]*phong_component*this->color[0]);
        temp_2=(light->color[1]*lambert_component*this->color[1])+(light->color[1]*phong_component*this->color[1]);
        temp_3=(light->color[2]*lambert_component*this->color[2])+(light->color[2]*phong_component*this->color[2]);
        color_pointer_quad[0]=color_pointer_quad[0]+temp_1;
        color_pointer_quad[1]=color_pointer_quad[1]+temp_2;
        color_pointer_quad[2]=color_pointer_quad[2]+temp_3;
    }
    void computeRecursiveColor(double* color_pointer_quad,double* ref_col )
    {
        for(int i=0; i<3; i++)
        {
            color_pointer_quad[i]=color_pointer_quad[i]+ref_col[i]*this->coEfficients[3];
        }
    }
    void computePhong_and_lambert(Ray inc,Vector3D normal,Ray *r, Light *light,double* color_pointer_sphere)
    {
        double lambert=dot(inc.Rd,normal);
        Vector3D reflected=reflected_ray(&inc,normal);
        double phong=dot(reflected,r->Rd);
        this->computeReflectionColor(color_pointer_sphere,light,lambert,phong);
    }
    Vector3D change_normal(Ray *r,Vector3D normal)
    {
        Vector3D temp=multiplyScalar(r->Rd,-1.0);
        //return normal;
        if(dot(temp,normal)>0.0) return normal;
        else
        {
            double x1=-normal.x;
            double y1=-normal.y;
            double z1=-normal.z;
            Vector3D norm_temp(x1,y1,z1);
            return norm_temp;
        }
        //else return multiplyScalar(normal,-1.0);
    }
    bool check_if_inside(Vector3D intersection_point)
    {
        if(length!=0 && (intersection_point.x<reference_point.x||intersection_point.x>reference_point.x+length))
        {
            return false;
        }
        if(width!=0 && (intersection_point.y<reference_point.y||intersection_point.y>reference_point.y+width))
        {
            return false;
        }
        if(height!=0 && (intersection_point.z<reference_point.z||intersection_point.z>reference_point.z+height))
        {
            return false;
        }
        return true;

    }
    double get_parameter_t(Ray* r)
    {
        double a,b,c,dis;
        a = parameters.get_a()*pow(r->Rd.x,2)+parameters.get_b()* pow(r->Rd.y,2)+ parameters.get_c()*pow(r->Rd.z,2)+parameters.get_d()*r->Rd.x*r->Rd.y
            + parameters.get_f()*r->Rd.y*r->Rd.z+parameters.get_e()*r->Rd.x *r->Rd.z;

        b=2.0*parameters.get_a()*r->Ro.x*r->Rd.x+2*parameters.get_b()*r->Ro.y*r->Rd.y+
          2.0*parameters.get_c()*r->Ro.z*r->Rd.z+parameters.get_d()*((r->Ro.x * r->Rd.y) + (r->Ro.y * r->Rd.x))+parameters.get_e()*(r->Ro.x * r->Rd.z + r->Ro.z * r->Rd.x)+parameters.get_f()*(r->Ro.y* r->Rd.z + r->Rd.y * r->Ro.z)
          + parameters.get_h()* r->Rd.y + parameters.get_i() * r->Rd.z+parameters.get_g()*r->Rd.x;

        c = parameters.get_a()*pow(r->Ro.x,2) + parameters.get_b()*pow(r->Ro.y,2)+parameters.get_c()*pow(r->Ro.z,2)+parameters.get_d()*r->Ro.x * r->Ro.y + parameters.get_e() * r->Ro.x * r->Ro.z + parameters.get_f() * r->Ro.y * r->Ro.z +
            parameters.get_g()*r->Ro.x +parameters.get_j()+parameters.get_h()*r->Ro.y + parameters.get_i() * r->Ro.z;
        if(a==0)
        {
            return ((-c/b)*(1.0));
        }

        dis = (b*b)-(4.0*a*c);
        if(dis<0) return -1;
        double root_1,root_2;
        double bval=-b;
        double nume=(2.0*a);
        root_1=(bval-sqrt(dis))/nume;


        //if (root_1<0.0) root_1=-9999;
        Vector3D intersection_root1((r->Ro.x+(r->Rd.x*root_1)),(r->Ro.y+(r->Rd.y*root_1)),(r->Ro.z+(r->Rd.z*root_1)));
        if(check_if_inside(intersection_root1)==false) root_1=-9999;


        root_2=(bval+sqrt(dis))/nume;
        Vector3D intersection_root2((r->Ro.x+(r->Rd.x*root_2)),(r->Ro.y+(r->Rd.y*root_2)),(r->Ro.z+(r->Rd.z*root_2)));
        if(check_if_inside(intersection_root2)==false) root_2=-9999;

        if(root_1!=-9999)
        {
            if(root_2!=-9999) return min(root_1,root_2);
            else return root_1;
        }
        else
        {
            if(root_2!=-9999) return root_2;
            else return -1;
        }

    }
    double intersect(Ray *,double *,int);

};
double General_quadratic::intersect(Ray *r, double *color_pointer_quad, int level)
{
    double t;
    t=get_parameter_t(r);
    if(t>0)
    {
        if(level==0)
        {
            return t;
        }
        Vector3D intersection_point((r->Ro.x+(r->Rd.x*t)),(r->Ro.y+(r->Rd.y*t)),(r->Ro.z+(r->Rd.z*t)));
        double norm_x,norm_y,norm_z;
        norm_x=2.0*parameters.get_a()*intersection_point.x+parameters.get_d()*intersection_point.y;
        norm_x=norm_x+parameters.get_e()*intersection_point.z+parameters.get_g();
        norm_y=2.0*parameters.get_b()*intersection_point.y+parameters.get_d()*intersection_point.x;
        norm_y=norm_y+parameters.get_f()*intersection_point.z+parameters.get_h();
        norm_z=2.0*parameters.get_c()*intersection_point.z+parameters.get_e()*intersection_point.x;
        norm_z=norm_z+parameters.get_f()*intersection_point.y+parameters.get_i();
        Vector3D normal(norm_x,norm_y,norm_z);
        computeAmbianceColor(color_pointer_quad);
        normal=normalize(normal);
        normal=change_normal(r,normal);

        for(int i=0; i<lights.size(); i++)
        {

            Vector3D light_dir=createNewPoint_Sub(lights[i]->light_position,intersection_point);
            Vector3D temp=light_dir;

            light_dir=normalize(light_dir);
            if(lights[i]->cut_off_angle!=INF)
            {
                Vector3D temp=multiplyScalar(light_dir,-1.0);
                double angle=check_angle(temp,lights[i]->direction);
                if(angle>lights[i]->cut_off_angle) continue;
            }
            Vector3D light_pos=get_light_pos(light_dir,intersection_point);
            Ray incident_ray(light_pos,light_dir);

            double length_of_light_dir=sqrt((temp.x*temp.x)+(temp.y*temp.y)+(temp.z*temp.z));
            int flag=check_if_obscured(incident_ray,length_of_light_dir);
            //check if ray is obscured
            if(flag==0)
            {
                 double lambert=dot(incident_ray.Rd,normal);
                 Vector3D reflected=reflected_ray(&incident_ray,normal);
                 double phong=dot(reflected,r->Rd);
                 computeReflectionColor(color_pointer_quad,lights[i],lambert,phong);
            }

        }

        if(level>=level_of_rec) return t;

        else
        {

            Vector3D refl_again=reflected_ray(r,normal);
            Vector3D refl_pos=get_light_pos(refl_again,intersection_point);
            Ray Ref_ray(refl_pos,refl_again);
            double tMin=INF;
            double t;
            int nearest=-9999;
            for(int k=0; k<objects.size(); k++)
            {
                double *ref_Color=new double[3];
                ref_Color[0]=0.0;
                ref_Color[1]=0.0;
                ref_Color[2]=0.0;
                t=objects[k]->intersect(&Ref_ray,ref_Color,0);
                if(check_t(t,tMin))
                {

                    nearest=k;
                    tMin=t;

                }
                delete_array(ref_Color);
            }
            if(nearest!=-9999)
            {
                double *ref_Color=new double[3];
                ref_Color[0]=0.0;
                ref_Color[1]=0.0;
                ref_Color[2]=0.0;
                tMin=objects[nearest]->intersect(&Ref_ray,ref_Color,level+1);
                if(tMin!=-1)
                {
                    computeRecursiveColor(color_pointer_quad,ref_Color);
                }
                delete_array(ref_Color);

            }
            color_pointer_quad[0]>1?color_pointer_quad[0]=1:(color_pointer_quad[0]<0?color_pointer_quad[0]=0:color_pointer_quad[0]=color_pointer_quad[0]);
            color_pointer_quad[1]>1?color_pointer_quad[1]=1:(color_pointer_quad[1]<0?color_pointer_quad[1]=0:color_pointer_quad[1]=color_pointer_quad[1]);
            color_pointer_quad[2]>1?color_pointer_quad[2]=1:(color_pointer_quad[2]<0?color_pointer_quad[2]=0:color_pointer_quad[2]=color_pointer_quad[2]);
            //cout<<tMin<<endl;
            return tMin;
        }
    }
    else return -1;
}

