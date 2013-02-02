#ifndef GLWINDOW_H
#define GLWINDOW_H

#include <GL/glut.h>

const int SCREEN_WIDTH;
const int SCREEN_HEIGHT;

void draw();
void myInit();
//void initMenu();
void update();
void resize(int, int);
void keyHandler(unsigned char, int, int);
void specialKey(int, int, int);
void computeAngles();
void menuHandler();

#endif
