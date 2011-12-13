#include <stdlib.h>
#include <stdio.h>
#include <GL/freeglut.h>

#define XSIZE 96
#define YSIZE 68
#define SCALE 4

typedef unsigned char u8;

static u8* pPixelBuf;
static u8* pDisplayBuf;

float posX = 65536.0;
float posY = 65536.0;

typedef float Vector[3];

float time = 0;

void setDispPixel(int _x, int _y, u8 value)
{
	int x,y;
	for (x = 0; x < SCALE; x++)
		for (y = 0; y < SCALE; y++)
		{
			int pos = (_y*SCALE+y)*XSIZE*SCALE+_x*SCALE+x;
			pos += 2*pos;
			int i = 0;
			for (i = 0; i < 3; i++)
			{
				pDisplayBuf[pos+i] = value;  
				if ((x % SCALE == 0) || (y % SCALE ==0))
					pDisplayBuf[pos+i] /= 2;
			}
		}
}

void setPixel(int x, int y, u8 value)
{
	pPixelBuf[y*XSIZE+x] = value;
}

void addPixel(int x, int y, u8 value)
{
	pPixelBuf[y*XSIZE+x] = value;
}

float sphereIntersection(Vector spherePos, Vector rayOrg, Vector rayDir)
{
	#define RADIUS 1.6

	float a,b,c,t = -100.0;

	float ocVec[3];
	ocVec[0] = rayOrg[0] - spherePos[0];
	ocVec[1] = rayOrg[1] - spherePos[1];
	ocVec[2] = rayOrg[2] - spherePos[2];

	a = rayDir[0]*rayDir[0] + rayDir[1]*rayDir[1] + rayDir[2]*rayDir[2];
	b = 2*(ocVec[0]*rayDir[0]+ocVec[1]*rayDir[1]+ocVec[2]*rayDir[2]);
	c = ocVec[0]*ocVec[0]+ocVec[1]*ocVec[1]+ocVec[2]*ocVec[2]  - (RADIUS*RADIUS);

	float disc = b * b - 4 * a * c;

#define BORDER 0.25f
	if (disc >= -BORDER && disc <= BORDER)
	{
		return 0;
	}

	if (disc < -BORDER) { //setPixel(x,y,0); 
		return t;

	}

	float distSqrt = sqrtf(disc);
	float q;
	if (b < 0)
		q = (-b - distSqrt)/2.0;
	else
		q = (-b + distSqrt)/2.0;

    // compute t0 and t1
    float t0 = q / a;
    float t1 = c / q;

    // make sure t0 is smaller than t1
    if (t0 > t1)
    {
        // if t0 is bigger than t1 swap them around
        float temp = t0;
        t0 = t1;
        t1 = temp;
    }

    // if t1 is less than zero, the object is in the ray's negative direction
    // and consequently the ray misses the sphere
   if (t1 < 0) return t1;


    // if t0 is less than zero, the intersection point is at t1
    if (t0 < 0)
    {
        t = t1;
    }
    // else the intersection point is at t0
    else
    {
        t = t0;
    }


//printf("%f \t%f \t%f \t --- %f \t%f \t%f \t --- %f\n",rayOrg[0],rayOrg[1],rayOrg[2],rayDir[0],rayDir[1],rayDir[2],t);
	return t;
}

float planeIntersection(Vector rayOrg, Vector rayDir)
{
	 return -(4.0-rayOrg[1])/(rayDir[1]);
}

u8 checkerBoard(float t, float x, float z)
{
	float dir = (t > 0) ? 1 : -1;
	float xC = (x*0.15) + dir*posX;
	float zC = (z*0.15) - dir*posY;

	int checker = ((int)xC % 2 + (int)zC%2) % 2;

	int dist = 128 - (int)sqrt(x*x + z*z)*2;
	if (dist < 0) dist = 0; 

	return dist*(int)checker;
}


void draw()
{
	int x,y;
	float rayOrg[3],rayDir[3];

	Vector spherePos;
	spherePos[1] = 0.0f;
	spherePos[2] = 2*cos(time/1600)+1;
	spherePos[0] = 1.2*sin(time/400);


	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
			setPixel(x,y,0);

	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
		{
			float t = 0;
			rayOrg[0] = 0;
			rayOrg[1] = 0;
			rayOrg[2] = -3.5;

			float rayDir[3];
			rayDir[0] = 1.5*((float)x/(float)(XSIZE)-0.5f);
			rayDir[1] = 1.5*((float)y/(float)(YSIZE)*0.5f*(float)XSIZE/(float)YSIZE-0.33f);
			rayDir[2] = 1;

			t = sphereIntersection(spherePos,rayOrg,rayDir);
			if (t < 0) 
			{
				t = planeIntersection(rayOrg,rayDir);
				rayOrg[0] += t*rayDir[0];
				rayOrg[2] += t*rayDir[2];
				setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
			} else/*
			if (t == 0)
			{
				setPixel(x,y,128);
			}
			else*/
			{
				rayOrg[0] = rayOrg[0] + t*rayDir[0];
				rayOrg[1] = rayOrg[1] + t*rayDir[1];
				rayOrg[2] = rayOrg[2] + t*rayDir[2];

				Vector ocVec;
				ocVec[0] = rayOrg[0] - spherePos[0];
				ocVec[1] = rayOrg[1] - spherePos[1];
				ocVec[2] = rayOrg[2] - spherePos[2];
		//		printf("%f \t%f \t%f \t --- %f \t%f \t%f \t --- %f\n",rayOrg[0],rayOrg[1],rayOrg[2],rayDir[0],rayDir[1],rayDir[2],t);

				float reflProd = 2*(ocVec[0]*rayDir[0] + ocVec[1]*rayDir[1] + ocVec[2]*rayDir[2]);
				rayDir[0] = rayDir[0] - reflProd * ocVec[0] ;
				rayDir[1] = rayDir[1] - reflProd * ocVec[1] ;
				rayDir[2] = rayDir[2] - reflProd * ocVec[2] ;

				t = planeIntersection(rayOrg,rayDir);
				rayOrg[0] += t*rayDir[0];
				rayOrg[2] += t*rayDir[2];
				setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
			}
		}
}

void dither()
{
	int x,y;
	for (y = 0; y < YSIZE-1; y++)
	{
		for (x = 0; x < XSIZE-1; x++)
		{
			int oldPixel = pPixelBuf[y*XSIZE+x];
			int newPixel = oldPixel/128*128;

			int quantError = oldPixel - newPixel;
		
			pPixelBuf[y*XSIZE+x+1    ] += quantError / 2;
			pPixelBuf[(y+1)*XSIZE+x-1] += quantError / 4;
			pPixelBuf[(y+1)*XSIZE+x  ] += quantError / 4;
		
			pPixelBuf[y*XSIZE+x] = newPixel/128*255;
		}
	}

}


void display( )
{
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);

	draw();
	dither();

	int x,y;
	for (y = 0; y < YSIZE-1; y++)
		for (x = 0; x < XSIZE-1; x++) 
			setDispPixel(x,y,pPixelBuf[y*XSIZE+x]);

	
	glDrawPixels(XSIZE*SCALE, YSIZE*SCALE, GL_RGB, GL_UNSIGNED_BYTE, pDisplayBuf);

	glutSwapBuffers();
}

void update(int value)
{
	glutPostRedisplay();

	float dX = 0.07*sin(time/2000), dY = 0.05;
	float l = sqrt(dX*dX + dY*dY); 
	dX /= l*12; dY /= l*12;

	posX += dX;
	posY += dY;


	time += 1000/24;

    glutTimerFunc(1000/24, update, 0);
}

void main(int ac, char* av[])
{
// Init OpenGL stuff ///////////////////////////////////
    glutInit(&ac, av);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(XSIZE*SCALE, YSIZE*SCALE);
    glutCreateWindow("W1l510n 0R30 --- r0ket++ ");
    glutDisplayFunc(&display);


  	glutTimerFunc(1000/24, update, 0);

	pPixelBuf = (u8*)malloc( XSIZE*YSIZE*3 );
	pDisplayBuf = (u8*)malloc( XSIZE*SCALE * YSIZE*SCALE * 3 );

	glutMainLoop();
	free(pDisplayBuf);

	return EXIT_SUCCESS;
}


