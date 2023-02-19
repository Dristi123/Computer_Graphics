#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
};

struct point pos;
struct point u;
struct point r;
struct point l;
double square_side=20.0;
double sphere_radius=15.0;
void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}

void moveForward(double unit) {
    pos.x=pos.x+unit*l.x;
    pos.y=pos.y+unit*l.y;
    pos.z=pos.z+unit*l.z;
    }
void moveBackward(double unit) {
    pos.x=pos.x-unit*l.x;
    pos.y=pos.y-unit*l.y;
    pos.z=pos.z-unit*l.z;
    }
void moveRight(double unit) {
    pos.x=pos.x+unit*r.x;
    pos.y=pos.y+unit*r.y;
    pos.z=pos.z+unit*r.z;
    }
void moveLeft(double unit) {
    pos.x=pos.x-unit*r.x;
    pos.y=pos.y-unit*r.y;
    pos.z=pos.z-unit*r.z;
    }
void moveUp(double unit) {
    pos.x=pos.x+unit*u.x;
    pos.y=pos.y+unit*u.y;
    pos.z=pos.z+unit*u.z;
    }
void moveDown(double unit) {
    pos.x=pos.x-unit*u.x;
    pos.y=pos.y-unit*u.y;
    pos.z=pos.z-unit*u.z;
    }
void lookLeft(double unit) {
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(u.y*l.z-u.z*l.y)*sin(angle);
    l.y=l.y*cos(angle)+(u.z*l.x-u.x*l.z)*sin(angle);
    l.z=l.z*cos(angle)+(u.x*l.y-u.y*l.x)*sin(angle);
    r.x=r.x*cos(angle)+(u.y*r.z-u.z*r.y)*sin(angle);
    r.y=r.y*cos(angle)+(u.z*r.x-u.x*r.z)*sin(angle);
    r.z=r.z*cos(angle)+(u.x*r.y-u.y*r.x)*sin(angle);

    }
void lookRight(double unit) {
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(u.y*l.z-u.z*l.y)*sin(-angle);
    l.y=l.y*cos(angle)+(u.z*l.x-u.x*l.z)*sin(-angle);
    l.z=l.z*cos(angle)+(u.x*l.y-u.y*l.x)*sin(-angle);
    r.x=r.x*cos(angle)+(u.y*r.z-u.z*r.y)*sin(-angle);
    r.y=r.y*cos(angle)+(u.z*r.x-u.x*r.z)*sin(-angle);
    r.z=r.z*cos(angle)+(u.x*r.y-u.y*r.x)*sin(-angle);
    }
void lookUp(double unit) {
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(r.y*l.z-r.z*l.y)*sin(angle);
    l.y=l.y*cos(angle)+(r.z*l.x-r.x*l.z)*sin(angle);
    l.z=l.z*cos(angle)+(r.x*l.y-r.y*l.x)*sin(angle);
    u.x=u.x*cos(angle)+(r.y*u.z-r.z*u.y)*sin(angle);
    u.y=u.y*cos(angle)+(r.z*u.x-r.x*u.z)*sin(angle);
    u.z=u.z*cos(angle)+(r.x*u.y-r.y*u.x)*sin(angle);
    }
void lookDown(double unit) {
    angle=(unit*pi/180.0);
    l.x=l.x*cos(angle)+(r.y*l.z-r.z*l.y)*sin(-angle);
    l.y=l.y*cos(angle)+(r.z*l.x-r.x*l.z)*sin(-angle);
    l.z=l.z*cos(angle)+(r.x*l.y-r.y*l.x)*sin(-angle);
    u.x=u.x*cos(angle)+(r.y*u.z-r.z*u.y)*sin(-angle);
    u.y=u.y*cos(angle)+(r.z*u.x-r.x*u.z)*sin(-angle);
    u.z=u.z*cos(angle)+(r.x*u.y-r.y*u.x)*sin(-angle);
    }
void tiltClock(double unit) {
    angle=(unit*pi/180.0);
    u.x=u.x*cos(angle)+(l.y*u.z-l.z*u.y)*sin(angle);
    u.y=u.y*cos(angle)+(l.z*u.x-l.x*u.z)*sin(angle);
    u.z=u.z*cos(angle)+(l.x*u.y-l.y*u.x)*sin(angle);
    r.x=r.x*cos(angle)+(l.y*r.z-l.z*r.y)*sin(angle);
    r.y=r.y*cos(angle)+(l.z*r.x-l.x*r.z)*sin(angle);
    r.z=r.z*cos(angle)+(l.x*r.y-l.y*r.x)*sin(angle);
    }
void tiltCounterclock(double unit) {
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
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    glColor3f(1,1,1);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,0);
		glVertex3f( a,-a,0);
		glVertex3f(-a,-a,0);
		glVertex3f(-a, a,0);
	}glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
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
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
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
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
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
			}glEnd();
		}
	}
}
void drawPartialphere(double radius,int slices,int stacks) {
    struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*0.50*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*0.50*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        glColor3f(1,0,0);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                //lower hemisphere
			}glEnd();
		}
	}
    }
void drawPartialcylinder(double radius,double height,int slices,int stacks) {
    struct point points[100];
    int i=0;
    for(i=0;i<=slices;i++) {
        points[i].x=radius*cos(((double)i/(double)slices)*0.5*pi);
        points[i].y=radius*sin(((double)i/(double)slices)*0.5*pi);
        points[i].z=height;
    }
    int j=0;
    for(j=0;j<slices;j++) {
        glColor3f(0,1,0);
        glBegin(GL_QUADS);

        glVertex3f(points[j].x, points[j].y, points[j].z);
        glVertex3f(points[j + 1].x, points[j + 1].y, points[j + 1].z);
        glVertex3f(points[j + 1].x, points[j + 1].y, -points[j + 1].z);
        glVertex3f(points[j].x, points[j].y, -points[j].z);

        glEnd();
    }

    }

void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

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


void specialKeyListener(int key, int x,int y){
	switch(key){
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
		    if(square_side>0){
		    square_side=square_side-1;
		    sphere_radius=sphere_radius+1;
		    }
			break;
		case GLUT_KEY_END:
		    if(sphere_radius>0){
		    square_side=square_side+1;
		    sphere_radius=sphere_radius-1;
		    }
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
				drawaxes=1-drawaxes;
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
void draw_spheres() {
    glPushMatrix();
    int i=0;
    glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    for(i=0;i<3;i++) {
    //glPopMatrix();
	glTranslated(-square_side,-square_side,-square_side);
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    }
    glPopMatrix();
    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
	for(i=0;i<3;i++) {
    //glPopMatrix();
	glTranslated(-square_side,-square_side,-square_side);
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    }
    glPopMatrix();
    }
void draw_spheres2() {
    glPushMatrix();

    glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(180,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(-90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(180,1,0,0);
    glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
	glPopMatrix();
	glPushMatrix();
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(180,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    glPushMatrix();
    glRotated(-90,0,0,1);
	glTranslated(square_side,square_side,square_side);
	drawPartialphere(sphere_radius,50,70);
    glPopMatrix();
    }
void draw_cylinders() {
    glPushMatrix();
    int i=0;
    glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    for(i=0;i<3;i++) {
    //glPopMatrix();
	glTranslated(-square_side,-square_side,0);
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    }
    glPopMatrix();
    glPushMatrix();
    glRotated(90,1,0,0);
    int j=0;
    glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    for(j=0;j<3;j++) {
    //glPopMatrix();
	glTranslated(-square_side,-square_side,0);
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    }
    glPopMatrix();
    glPushMatrix();
    glRotated(90,0,1,0);
    int k=0;
    glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    for(k=0;k<3;k++) {
    //glPopMatrix();
	glTranslated(-square_side,-square_side,0);
    glRotated(90,0,0,1);
	glTranslated(square_side,square_side,0);
	drawPartialcylinder(sphere_radius,square_side,50,70);
    }
    glPopMatrix();
    }
void draw_squares() {
    glPushMatrix();
    glTranslated(0,0,square_side+sphere_radius);
    drawSquare(square_side);
    glPopMatrix();
    glPushMatrix();
    glTranslated(0,0,-(square_side+sphere_radius));
    drawSquare(square_side);
    glPopMatrix();

    int i=0;
    for(i=0;i<2;i++) {
    glPushMatrix();
    if(i==0)glRotated(90,1,0,0);
    else glRotated(-90,1,0,0);
    glTranslated(0,0,square_side+sphere_radius);
    drawSquare(square_side);
    glPopMatrix();
    }
    int j=0;
    for(j=0;j<2;j++){
    glPushMatrix();
    if(j==0)glRotated(90,0,1,0);
    else glRotated(-90,0,1,0);
    glTranslated(0,0,square_side+sphere_radius);
    drawSquare(square_side);
    glPopMatrix();
    }
}
void draw_intermediate() {
    draw_spheres();
    draw_cylinders();
    draw_squares();
}


void display(){

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
	drawGrid();


    //drawSquare(10);

    //drawSS();
     draw_intermediate();
    //drawCircle(30,24);

    //drawCone(20,50,24);
    //glColor3f(1,1,0);






	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=0.05;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;
    pos.x=100.0;
    pos.y=100.0;
    pos.z=0.0;
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
	gluPerspective(80,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	init();

	glEnable(GL_DEPTH_TEST);	//enable Depth Testing

	glutDisplayFunc(display);	//display callback function
	glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

	glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);
	glutMouseFunc(mouseListener);

	glutMainLoop();		//The main loop of OpenGL

	return 0;
}
