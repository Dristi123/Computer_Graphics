#include <windows.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "1705029_Header.h"
#include "bitmap_image.hpp"
#include<bits/stdc++.h>
#include<stdio.h>
#define INF numeric_limits<double>::infinity()
#include<math.h>
#include<stdlib.h>
#include<iostream>
#define pi (2*acos(0.0))
using namespace std;

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle,angle2,angle3,angle4;
int bitmap_count=0;
int fov=80;


extern vector<Object*> objects;
extern vector<Light*> lights;
extern int level_of_rec;
int image_dim;


extern struct Vector3D pos;
Vector3D u;
Vector3D r;
Vector3D l;
int window_height=500;
int window_width=500;

void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 0.0, 0.0);
        glBegin(GL_LINES);
        {
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);
            glColor3f(0.0, 1.0, 0.0);
            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);
            glColor3f(0.0, 0.0, 1.0);
            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }
        glEnd();
    }
}
void moveForward(double unit)
{
    pos.x=pos.x+unit*l.x;
    pos.y=pos.y+unit*l.y;
    pos.z=pos.z+unit*l.z;
}

void moveBackward(double unit)
{
    pos.x=pos.x-unit*l.x;
    pos.y=pos.y-unit*l.y;
    pos.z=pos.z-unit*l.z;
}
void moveRight(double unit)
{
    pos.x=pos.x+unit*r.x;
    pos.y=pos.y+unit*r.y;
    pos.z=pos.z+unit*r.z;
}
void moveLeft(double unit)
{
    pos.x=pos.x-unit*r.x;
    pos.y=pos.y-unit*r.y;
    pos.z=pos.z-unit*r.z;
}
void moveUp(double unit)
{
    pos.x=pos.x+unit*u.x;
    pos.y=pos.y+unit*u.y;
    pos.z=pos.z+unit*u.z;
}
void moveDown(double unit)
{
    pos.x=pos.x-unit*u.x;
    pos.y=pos.y-unit*u.y;
    pos.z=pos.z-unit*u.z;
}
void lookLeft(double unit)
{
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(u.y*l.z-u.z*l.y)*sin(angle);
    l.y=l.y*cos(angle)+(u.z*l.x-u.x*l.z)*sin(angle);
    l.z=l.z*cos(angle)+(u.x*l.y-u.y*l.x)*sin(angle);
    r.x=r.x*cos(angle)+(u.y*r.z-u.z*r.y)*sin(angle);
    r.y=r.y*cos(angle)+(u.z*r.x-u.x*r.z)*sin(angle);
    r.z=r.z*cos(angle)+(u.x*r.y-u.y*r.x)*sin(angle);

}
void lookRight(double unit)
{
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(u.y*l.z-u.z*l.y)*sin(-angle);
    l.y=l.y*cos(angle)+(u.z*l.x-u.x*l.z)*sin(-angle);
    l.z=l.z*cos(angle)+(u.x*l.y-u.y*l.x)*sin(-angle);
    r.x=r.x*cos(angle)+(u.y*r.z-u.z*r.y)*sin(-angle);
    r.y=r.y*cos(angle)+(u.z*r.x-u.x*r.z)*sin(-angle);
    r.z=r.z*cos(angle)+(u.x*r.y-u.y*r.x)*sin(-angle);
}
void lookUp(double unit)
{
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(r.y*l.z-r.z*l.y)*sin(angle);
    l.y=l.y*cos(angle)+(r.z*l.x-r.x*l.z)*sin(angle);
    l.z=l.z*cos(angle)+(r.x*l.y-r.y*l.x)*sin(angle);
    u.x=u.x*cos(angle)+(r.y*u.z-r.z*u.y)*sin(angle);
    u.y=u.y*cos(angle)+(r.z*u.x-r.x*u.z)*sin(angle);
    u.z=u.z*cos(angle)+(r.x*u.y-r.y*u.x)*sin(angle);
}
void lookDown(double unit)
{
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(r.y*l.z-r.z*l.y)*sin(-angle);
    l.y=l.y*cos(angle)+(r.z*l.x-r.x*l.z)*sin(-angle);
    l.z=l.z*cos(angle)+(r.x*l.y-r.y*l.x)*sin(-angle);
    u.x=u.x*cos(angle)+(r.y*u.z-r.z*u.y)*sin(-angle);
    u.y=u.y*cos(angle)+(r.z*u.x-r.x*u.z)*sin(-angle);
    u.z=u.z*cos(angle)+(r.x*u.y-r.y*u.x)*sin(-angle);
}
void tiltClock(double unit)
{
    angle=(unit*pi/180.0);
    u.x=u.x*cos(angle)+(l.y*u.z-l.z*u.y)*sin(angle);
    u.y=u.y*cos(angle)+(l.z*u.x-l.x*u.z)*sin(angle);
    u.z=u.z*cos(angle)+(l.x*u.y-l.y*u.x)*sin(angle);
    r.x=r.x*cos(angle)+(l.y*r.z-l.z*r.y)*sin(angle);
    r.y=r.y*cos(angle)+(l.z*r.x-l.x*r.z)*sin(angle);
    r.z=r.z*cos(angle)+(l.x*r.y-l.y*r.x)*sin(angle);
}
void tiltCounterclock(double unit)
{
    angle=(unit*pi/180.0);
    u.x=u.x*cos(angle)+(l.y*u.z-l.z*u.y)*sin(-angle);
    u.y=u.y*cos(angle)+(l.z*u.x-l.x*u.z)*sin(-angle);
    u.z=u.z*cos(angle)+(l.x*u.y-l.y*u.x)*sin(-angle);
    r.x=r.x*cos(angle)+(l.y*r.z-l.z*r.y)*sin(-angle);
    r.y=r.y*cos(angle)+(l.z*r.x-l.x*r.z)*sin(-angle);
    r.z=r.z*cos(angle)+(l.x*r.y-l.y*r.x)*sin(-angle);
}

void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-8; i<=8; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);

                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);
            }
        }
        glEnd();
    }
}
void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    glColor3f(1,1,1);
    glBegin(GL_QUADS);
    {
        glVertex3f( a, a,0);
        glVertex3f( a,-a,0);
        glVertex3f(-a,-a,0);
        glVertex3f(-a, a,0);
    }
    glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct Vector3D points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0; i<=segments; i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0; i<segments; i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct Vector3D points[100];
    //generate points
    for(i=0; i<=segments; i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0; i<segments; i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
    struct Vector3D points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0; i<=stacks; i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
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
        glColor3f(1,0,0);
        for(j=0; j<slices; j++)
        {
            glBegin(GL_QUADS);
            {
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }
            glEnd();
        }
    }
}


bool check_tval(double t,double min_t)
{
    if(t>0 && t<min_t) return true;
    else return false;
}
void specialKeyListener(int key, int x,int y)
{
    switch(key)
    {
    case GLUT_KEY_DOWN:		//down arrow key
        cameraHeight -= 3.0;
        moveBackward(3.0);
        break;
    case GLUT_KEY_UP:		// up arrow key
        cameraHeight += 3.0;
        moveForward(3.0);
        break;

    case GLUT_KEY_RIGHT:
        cameraAngle += 0.03;
        moveRight(3.0);
        break;
    case GLUT_KEY_LEFT:
        cameraAngle -= 0.03;
        moveLeft(3.0);
        break;

    case GLUT_KEY_PAGE_UP:
        moveUp(3.0);
        break;
    case GLUT_KEY_PAGE_DOWN:
        moveDown(3.0);
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:

        break;
    case GLUT_KEY_END:

        break;

    default:
        break;
    }
}

void capture()
{
    cout<<"Please wait.bitmap image is being captured"<<endl;
    ofstream check_file;
    check_file.open("S:\\STUDY\\CSE 410\\ofl9 3\\dristi\\check_mine.txt");
    bitmap_image image(image_dim,image_dim);

    int nearest=-9999;
    double du,dv,t;
    du=(double(window_width)/image_dim);
    dv=(double(window_height)/image_dim);
    double x_temp,y_temp,z_temp;
    double plane_distance;
    plane_distance=(window_height*0.5)/(tan(fov*0.5*(pi/180.0)));
    x_temp = pos.x +(l.x*plane_distance)-r.x*(window_width*0.5);
    x_temp=x_temp+u.x*(window_height*0.5);
    x_temp = x_temp+(r.x*(du/2.0))-(u.x*(dv/2.0));
    y_temp = pos.y +(l.y * plane_distance)-r.y*(window_width*0.5);
    y_temp=y_temp+u.y*(window_height*0.5);
    y_temp = y_temp+(r.y*(du/2.0))-(u.y*(dv/2.0));
    z_temp = pos.z +(l.z * plane_distance)-r.z*(window_width*0.5);
    z_temp=z_temp+u.z*(window_height*0.5);
    z_temp = z_temp+(r.z*(du/2.0))-(u.z*(dv/2.0));
    Vector3D leftmost(x_temp,y_temp,z_temp);
    for(int i=0; i<image_dim; i++)
    {
        for(int j=0; j<image_dim; j++)
        {
            image.set_pixel(i,j,0,0,0);
        }
    }
    for(int i=0; i<image_dim; i++)
    {
        for(int j=0; j<image_dim; j++)
        {
            double ray_x,ray_y,ray_z;
            Vector3D current_Pixel((leftmost.x+(r.x*(i*du))-(u.x*j*dv)),(leftmost.y+(r.y*(i*du))-(u.y*j*dv)),(leftmost.z+(r.z*(i*du))-(u.z*j*dv)));
            ray_x=current_Pixel.x-pos.x;
            ray_y=current_Pixel.y-pos.y;
            ray_z=current_Pixel.z-pos.z;
            double t_min=INF;
            Ray new_ray(pos,normalize(Vector3D(ray_x,ray_y,ray_z)));


            for(int k=0; k<objects.size(); k++)
            {
                double *color=new double[3];
                color[0]=0.0;
                color[1]=0.0;
                color[2]=0.0;
                t=objects[k]->intersect(&new_ray,color,0);

                if(check_tval(t,t_min))
                {

                    nearest=k;
                    t_min=t;

                }
            }
            if(nearest!=-9999)
            {
                double *color=new double[3];
                color[0]=0.0;
                color[1]=0.0;
                color[2]=0.0;
                t_min=objects[nearest]->intersect(&new_ray,color,1);
                double pixel_color[3];
                for(int k=0; k<3; k++)
                {
                    pixel_color[k]=color[k]*255;
                }
                image.set_pixel(i,j,(pixel_color[0]),(pixel_color[1]),(pixel_color[2]));
                check_file<<color[0]*255<<" "<<color[1]*255<<" "<<color[2]*255<<endl;

            }

        }

    }
    bitmap_count=bitmap_count+1;
    string number;
    std::stringstream ss;
    ss<<bitmap_count;
    ss>>number;
    string output_path="S:\\STUDY\\CSE 410\\ofl9 3\\dristi\\Output_1"+number+".bmp";
    image.save_image(output_path);;
    cout<<"bitmap image captured successfully"<<endl;
}


void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {
    case '0':
        capture();
        break;
    case '1':
        lookLeft(3);
        break;
    case '2':
        lookRight(3);
        break;
    case '3':
        lookUp(3);
        break;
    case '4':
        lookDown(3);
        break;
    case '5':
        tiltClock(3);
        break;
    case '6':
        tiltCounterclock(3);
        break;

    default:
        break;
    }
}




void mouseListener(int button, int state, int x, int y) 	//x, y is the x-y of the screen (2D)
{
    switch(button)
    {
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN) 		// 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            //drawaxes=1-drawaxes;
        }
        break;

    case GLUT_RIGHT_BUTTON:
        //........
        break;

    case GLUT_MIDDLE_BUTTON:
        //........
        break;

    default:
        break;
    }
}

void drawSS()
{
    for(int i=0; i<objects.size(); i++)
    {
        objects[i]->draw();
    }
    for(int i=0; i<lights.size(); i++)
    {
        lights[i]->draw();
    }
}

void display()
{

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x,pos.y+l.y,pos.z+l.z,	u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    //drawGrid();


    //drawSphere(10,50,50);

    drawSS();
    //draw_intermediate();
    //drawCircle(30,24);
    //cout<<objects.size()<<endl;


    //drawCone(20,50,24);
    //glColor3f(1,1,0);






    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{


    //codes for any changes in Models, Camera
    angle+=0.05;
    glutPostRedisplay();
}

void init()
{
    //codes for initialization
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;
    Vector3D temp(130,110,40);
    pos=temp;
    u.x=0.0;
    u.y=0.0;
    u.z=1.0;
    r.x=-(1/sqrt(2));
    r.y=(1/sqrt(2));
    r.z=0.0;
    l.x=-(1/sqrt(2));
    l.y=-(1/sqrt(2));
    l.z=0.0;
    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(fov,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}
void loadData()
{
    int no_of_obj,no_of_point_light,no_of_spot_light;
    ifstream filein("S:\\STUDY\\CSE 410\\ofl9 3\\dristi\\scene.txt");
    if (!filein)
    {
        cout<<"Error in opening scene.txt"<<endl;
        exit(1);

    }
    else

    {
        int line_count=1;
        string lineOfFile;
        while (!filein.eof())
        {
            if(line_count==1)
            {
                getline(filein,lineOfFile);
                line_count++;

                std::istringstream iss(lineOfFile);
                iss>>level_of_rec;

            }
            if(line_count==2)
            {
                getline(filein,lineOfFile);
                line_count++;
                std::istringstream iss(lineOfFile);
                iss>>image_dim;

                getline(filein,lineOfFile);
                line_count++;

            }
            if(line_count==4)
            {
                getline(filein,lineOfFile);
                line_count++;
                std::istringstream iss(lineOfFile);
                iss>>no_of_obj;
            }
            else
            {
                for(int j=0; j<no_of_obj; j++)
                {
                    getline(filein,lineOfFile);
                    line_count++;
                    cout<<lineOfFile<<endl;
                    if(lineOfFile.compare("sphere")==0)
                    {
                        cout<<"heree"<<endl;
                        double x,y,z;
                        int radius;
                        double r,g,b;
                        double amb,diff,specular,rec_ref_co;
                        int shine;
                        getline(filein,lineOfFile);
                        cout<<"ekhane lof"<<lineOfFile<<endl;
                        line_count++;
                        std::istringstream iss(lineOfFile);
                        iss>>x>>y>>z;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss2(lineOfFile);
                        iss2>>radius;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss3(lineOfFile);
                        iss3>>r>>g>>b;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss4(lineOfFile);
                        iss4>>amb>>diff>>specular>>rec_ref_co;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss5(lineOfFile);
                        iss5>>shine;
                        getline(filein,lineOfFile);
                        line_count++;
                        struct Vector3D center(x,y,z);
                        Object *sphere;
                        sphere=new Sphere(center,radius);
                        sphere->setColor(r,g,b);
                        sphere->setCoEfficients(amb,diff,specular,rec_ref_co);
                        sphere->setShine(shine);
                        objects.push_back(sphere);


                    }
                    else if(lineOfFile.compare("triangle")==0)
                    {
                        cout<<"ttttriangle"<<endl;
                        double x1,y1,z1,x2,y2,z2,x3,y3,z3;
                        double r,g,b;
                        double amb,diff,specular,rec_ref_co;
                        int shine;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss(lineOfFile);
                        iss>>x1>>y1>>z1;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss2(lineOfFile);
                        iss2>>x2>>y2>>z2;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss3(lineOfFile);
                        iss3>>x3>>y3>>z3;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss4(lineOfFile);
                        iss4>>r>>g>>b;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss5(lineOfFile);
                        iss5>>amb>>diff>>specular>>rec_ref_co;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss6(lineOfFile);
                        iss6>>shine;
                        getline(filein,lineOfFile);
                        line_count++;
                        Object *triangle;
                        triangle=new Triangle(x1,y1,z1,x2,y2,z2,x3,y3,z3);
                        triangle->setColor(r,g,b);
                        triangle->setCoEfficients(amb,diff,specular,rec_ref_co);
                        triangle->setShine(shine);
                        objects.push_back(triangle);

                    }
                    else if(lineOfFile.compare("general")==0)
                    {
                        cout<<"gennn"<<endl;
                        double a,b,c,d,e,f,g,h,i,j,x,y,z;
                        double length,height,width;
                        double amb,diff,specular,rec_ref_co;
                        double r,gr,bl;
                        int shine;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss(lineOfFile);
                        iss>>a>>b>>c>>d>>e>>f>>g>>h>>i>>j;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss2(lineOfFile);
                        iss2>>x>>y>>z>>length>>width>>height;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss4(lineOfFile);
                        iss4>>r>>gr>>bl;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss5(lineOfFile);
                        iss5>>amb>>diff>>specular>>rec_ref_co;
                        getline(filein,lineOfFile);
                        line_count++;
                        std::istringstream iss6(lineOfFile);
                        iss6>>shine;
                        getline(filein,lineOfFile);
                        line_count++;
                        //Vector3D ref_point(x,y,z);
                        GenCoeffs params(a,b,c,d,e,f,g,h,i,j);
                        Object *general;
                        general=new General_quadratic(params,x,y,z,length,width,height);
                        general->setColor(r,gr,bl);
                        general->setCoEfficients(amb,diff,specular,rec_ref_co);
                        general->setShine(shine);
                        objects.push_back(general);

                    }
                }
                getline(filein,lineOfFile);
                line_count++;
                std::istringstream iss_1(lineOfFile);
                iss_1>>no_of_point_light;
                for(int i=0; i<no_of_point_light; i++)
                {
                    double x,y,z,r,g,b;
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_1(lineOfFile);
                    iss_1>>x>>y>>z;
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_2(lineOfFile);
                    iss_2>>r>>g>>b;
                    //struct Vector3D pos(x,y,z);
                    Light *point_light;
                    point_light=new Light(x,y,z);
                    point_light->setColor(r,g,b);
                    point_light->setAngle(INF);
                    lights.push_back(point_light);


                }
                getline(filein,lineOfFile);
                line_count++;
                //cout<<"eta oi line"<<lineOfFile<<endl;
                getline(filein,lineOfFile);
                line_count++;
                std::istringstream iss_2(lineOfFile);
                iss_2>>no_of_spot_light;
                for(int i=0;i<no_of_spot_light;i++) {
                    double x,y,z,r,g,b;
                    double dir_x,dir_y,dir_z;
                    double cut_off;
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_1(lineOfFile);
                    iss_1>>x>>y>>z;
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_2(lineOfFile);
                    iss_2>>r>>g>>b;
                    //Vector3D pos(x,y,z);
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_3(lineOfFile);
                    iss_3>>dir_x>>dir_y>>dir_z;
                    Vector3D dir(dir_x,dir_y,dir_z);
                    getline(filein,lineOfFile);
                    line_count++;
                    std::istringstream iss_4(lineOfFile);
                    iss_4>>cut_off;
                    Light *spot;
                    spot=new Light(x,y,z);
                    spot->setColor(r,g,b);
                    spot->setDirection(dir);
                    printPoint(dir);
                    spot->setAngle(cut_off);
                    cout<<cut_off<<endl;
                    lights.push_back(spot);

                }
                break;

            }
        }
        cout<<level_of_rec<<endl;
        cout<<image_dim<<endl;
        cout<<no_of_obj<<endl;
        cout<<objects.size()<<endl;
        cout<<lights.size()<<endl;
    }


}

int main(int argc, char **argv)
{




    glutInit(&argc,argv);
    glutInitWindowSize(500,500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();
    loadData();
    Object *floor;
    floor=new Floor(1000,20);
    floor->setShine(17);
    floor->setCoEfficients(0.5,0.4,0.4,0.5);

    objects.push_back(floor);
    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();
    vector<Object*>().swap(objects);
    vector<Light*>().swap(lights);


    //The main loop of OpenGL

    return 0;
}
