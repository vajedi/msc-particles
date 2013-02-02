#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "ScalarCorr.h"
#include "anorm.h"
#include "simulate.h"

const int nTerms = 20 + 1;
const float gl_POINT_SIZE = 1.0f;
const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;

void myInit();
//void initMenu();
void update();
void resize(int, int);
void keyHandler(unsigned char, int, int);
void specialKey(int, int, int);
//void menuHandler();


int main(int argc, char* argv[])
{
	randomSeed();
	
	int mid = (nTerms - 1)/2;
	
	init(nTerms, mid);

	initParticles(nTerms, mid, SIGMA);
	//simulateParticles(a, nTerms, mid, SIGMA, 0.1);
	//testANorm(a, nTerms, mid);
	//testScalarCorr(a, nTerms, mid, SIGMA);
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);
    glutInitWindowPosition(10, 10);
    glutCreateWindow("Particles 2D");

    glutDisplayFunc(draw);
    glutKeyboardFunc(keyHandler);
    glutSpecialFunc(specialKey);
    glutReshapeFunc(resize);
    glutTimerFunc(15, update, 0);
/*    glutCreateMenu(menuHandler);
    glutAddMenuEntry("Full light", FullLight);
    glutAddMenuEntry("Darkness", Darkness);
    glutAddMenuEntry("Quit", Quit);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
*/    
    glPointSize( gl_POINT_SIZE );
    myInit();
    glutMainLoop();
	
    
    freeA( nTerms );
    return 0;
}

void myInit()
{
    glClearColor(1, 1, 1, 1);
    glEnable(GL_DEPTH_TEST);   // So that obscured objects aren't drawn.
    
    // set up projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    // 2D orthographic projection
    // x_min, x_max, y_min, y_max
    gluOrtho2D(0, 2, -1, 1);

    // set the modelview matrix to I
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glLineWidth(3);
    glPointSize(3.0f);
}

void resize(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(0.0, (double)w / (double)h, -2, 200.0);
    glMatrixMode(GL_MODELVIEW);
}

void keyHandler(unsigned char key, int x, int y)
{
    if (key == 'q' || key == 27)
        exit(0);
    glutPostRedisplay();
}

void specialKey(int key, int x, int y)
{
    switch(key)
    {
        case GLUT_KEY_UP:
            break;
        case GLUT_KEY_DOWN:
			break;
        case GLUT_KEY_LEFT:
            break;
        case GLUT_KEY_RIGHT:
            break;
    }
    glutPostRedisplay();
}
