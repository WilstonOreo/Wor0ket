// This file is part of wor0ket. wor0ket is free software: you can
// redistribute it and/or modify it under the terms of the GNU General Public
// License as published by the Free Software Foundation, version 2.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, write to the Free Software Foundation, Inc., 51
// Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// Copyright (2011), Wilston Oreo

#include <stdlib.h>
#include <stdio.h>
#include <GL/freeglut.h>

#define XSIZE 96 
#define YSIZE 68
#define SCALE 4

typedef unsigned char u8;

static u8* pPixelBuf;
static u8* pDisplayBuf;

int posX = 256;
int posY = 256;
int time = 0;

typedef int Veci3[3];

///////// Look up tables for required for integer sin() and integer cos()

static const u8 sin_table[] = {
0, 2, 4, 5, 7, 8, 10, 11, 13, 15, 16, 18, 19, 21, 22, 24, 
25, 27, 29, 30, 32, 33, 35, 36, 38, 39, 41, 43, 44, 46, 47, 49, 
50, 52, 53, 55, 56, 58, 59, 61, 62, 64, 65, 67, 69, 70, 72, 73, 
75, 76, 78, 79, 80, 82, 83, 85, 86, 88, 89, 91, 92, 94, 95, 97, 
98, 100, 101, 102, 104, 105, 107, 108, 110, 111, 112, 114, 115, 117, 118, 119, 
121, 122, 123, 125, 126, 128, 129, 130, 132, 133, 134, 136, 137, 138, 140, 141, 
142, 143, 145, 146, 147, 149, 150, 151, 152, 154, 155, 156, 157, 159, 160, 161, 
162, 163, 165, 166, 167, 168, 169, 171, 172, 173, 174, 175, 176, 177, 179, 180, 
181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 192, 193, 194, 195, 196, 197, 
198, 199, 200, 201, 202, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 
213, 213, 214, 215, 216, 217, 218, 218, 219, 220, 221, 222, 222, 223, 224, 225, 
225, 226, 227, 228, 228, 229, 230, 230, 231, 232, 232, 233, 234, 234, 235, 235, 
236, 237, 237, 238, 238, 239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 
245, 245, 245, 246, 246, 247, 247, 247, 248, 248, 249, 249, 249, 250, 250, 250, 
251, 251, 251, 251, 252, 252, 252, 253, 253, 253, 253, 253, 254, 254, 254, 254, 
254, 254, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

int intsin(int x)
{
	int sn = (x & 256) ? sin_table[255-(x & 255)] : sin_table[x & 255] ;
	if (x & 512) sn = -sn;
	return sn;
}

int intcos(int x)
{
	x += 256;
	int sn = (x & 256) ? sin_table[255-(x & 255)] : sin_table[x & 255] ;
	if (x & 512) sn = -sn;
	return sn;
}

/* by Mark Crowne 
 * http://www.azillionmonkeys.com/qed/sqroot.html
 * */
static unsigned int intsqrt (unsigned long val) {
  unsigned int temp, g=0;

	  if (val >= 0x40000000) {
	    g = 0x8000; 
	    val -= 0x40000000;
	  }

	#define INNER_ISQRT(s)                        \
	  temp = (g << (s)) + (1 << ((s) * 2 - 2));   \
	  if (val >= temp) {                          \
	    g += 1 << ((s)-1);                        \
	    val -= temp;                              \
	  }

	  INNER_ISQRT (15)
	  INNER_ISQRT (14)
	  INNER_ISQRT (13)
	  INNER_ISQRT (12)
	  INNER_ISQRT (11)
	  INNER_ISQRT (10)
	  INNER_ISQRT ( 9)
	  INNER_ISQRT ( 8)
	  INNER_ISQRT ( 7)
	  INNER_ISQRT ( 6)
	  INNER_ISQRT ( 5)
	  INNER_ISQRT ( 4)
	  INNER_ISQRT ( 3)
	  INNER_ISQRT ( 2)

	#undef INNER_ISQRT

	  temp = g+g+1;
	  if (val >= temp) g++;
	  return g;
}


///////// Set pixel functions

// Sets a pixel on the display buffer, with scale
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
				if (SCALE >= 4)
					if ((x % SCALE == 0) || (y % SCALE ==0))
						pDisplayBuf[pos+i] /= 2;
			}
			pDisplayBuf[pos+0] /= 2;
			pDisplayBuf[pos+2] /= 4;
		}
}

void setPixel(int x, int y, u8 value)
{
	pPixelBuf[y*XSIZE+x] = value;
}


// Ray-sphere intersection
int sphereIntersection(Veci3 spherePos, Veci3 rayOrg, Veci3 rayDir, int R)
{
	int t = -1000; 
	int ocVec[3];

	Veci3 rayDiri;
	int i;
	for (i = 0; i < 3; i++)
	{
		ocVec[i] = (rayOrg[i] - spherePos[i]) >> 8;
		rayDiri[i] = rayDir[i] >> 8;
	}

	// The equation to solve has the form
	// At² + Bt + C
	int A = rayDiri[0]*rayDiri[0] + rayDiri[1]*rayDiri[1] + rayDiri[2]*rayDiri[2];
	int B = ocVec[0]*rayDiri[0]+ocVec[1]*rayDiri[1]+ocVec[2]*rayDiri[2];
	int C = ocVec[0]*ocVec[0]+ocVec[1]*ocVec[1]+ocVec[2]*ocVec[2];

	A >>= 5; B >>= 4; C >>= 3;
	int disc = B * B - A* (C - ((R*R) >> 3));

	if (disc < 0)  return t;

	int distSqrt = intsqrt(disc);
	int q;
	if (B < 0)
		q = (-B - distSqrt) / 2;
	else
		q = (-B + distSqrt) / 2;

	int t0 = (!A) ? 1 << 30 : (q << 12) / A;
	int t1 = (!q) ? 1 << 30 : (C << 10) / q;

	if (t0 > t1)
	{
		int temp = t0;
		t0 = t1;
		t1 = temp;
	}

	// if t1 is less than zero, the object is in the ray's negative direction
	// and consequently the ray misses the sphere
	if (t1 < 0) return t1;

	// if t0 is less than zero, the intersection point is at t1
	if (t0 < 0) t = t1;
	else 		t = t0;

	return t;
}



int planeIntersection(Veci3 rayOrg, Veci3 rayDir)
{
	return (abs(rayDir[1]) < 4) ? 1 << 30 : -(((65536*4-rayOrg[1]) << 13)/rayDir[1]) >> 1;
}


// Draw a checker board, position changes over time
u8 checkerBoard(int t, int x, int z)
{
	int dir = (t > 0) ? 1 : -1;

	if (abs(z) > 65536*128) return 0;

	int pX = posX,pY = posY;
	int xC = ((x*96 >> 8) + dir*pX) >> 16;
	int zC = ((z*96 >> 8) - dir*pY) >> 16;
	int checker = ((u8)xC & 1) ^ ((u8)zC & 1);

	int dist = intsqrt((x >> 12)*(x >> 12)+ (z >> 12)*(z >> 12))/4-32;
	if (dist > 128) dist = 128;
	if (dist < 0) dist = 0;

	return (128-dist)*checker;
}


////////////////////////////////////
// Scene 1: Floor and ceiling with a checker board texture
// 	 		and a sphere in between reflecting them
void drawScene1()
{
	int x,y;
	Veci3 rayOrg, spherePos;

	spherePos[1] = 0;
	spherePos[2] = 256*2*intcos(time/8)+1;
	spherePos[0] = 256*intsin(time/2);

	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
		{
			int t = 0;
			rayOrg[0] = 0;
			rayOrg[1] = 0;
			rayOrg[2] = -65536*9/2;

			Veci3 rayDir;
			rayDir[0] = 2*((x << 16)/XSIZE - 32768); 
			rayDir[1] = 2*((y << 16)/YSIZE*XSIZE/YSIZE/2 - 24000   );
			rayDir[2] = 65536;
			t = sphereIntersection(spherePos,rayOrg,rayDir,512);

			if (t < 0) 
			{
				t = planeIntersection(rayOrg,rayDir);
				rayOrg[0] += (t >> 8)*(rayDir[0]) >> 4;
				rayOrg[2] += (t >> 8)*(rayDir[2]) >> 4;
				setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
			} else
				{
					rayOrg[0] = rayOrg[0]+ (t*rayDir[0] >> 12);
					rayOrg[1] = rayOrg[1]+ (t*rayDir[1] >> 12);
					rayOrg[2] = rayOrg[2]+ (t*rayDir[2] >> 12);

					Veci3 N;
					N[0] = (rayOrg[0] - spherePos[0]) >> 8;
					N[1] = (rayOrg[1] - spherePos[1]) >> 8;
					N[2] = (rayOrg[2] - spherePos[2]) >> 8;

					int reflProd = (N[0]*rayDir[0] + N[1]*rayDir[1] + N[2]*rayDir[2]) >> 9;
					rayDir[0] = rayDir[0] - (reflProd * N[0] >> 8);
					rayDir[1] = rayDir[1] - (reflProd * N[1] >> 8);
					rayDir[2] = rayDir[2] - (reflProd * N[2] >> 8);

					t = planeIntersection(rayOrg,rayDir);
					rayOrg[0] += (t >> 8)*(rayDir[0]) >> 4;
					rayOrg[2] += (t >> 8)*(rayDir[2]) >> 4;
					setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
				}
		}
}

////////////////////////////////////
// Scene 2: A moving sphere casting a shadow on the ground
//
void drawScene2()
{
	int x,y;
	Veci3 rayOrg, spherePos;

	spherePos[1] = 256*intsin(time/3);
	spherePos[2] = -512*intsin(time/5);
	spherePos[0] = 0;
	Veci3 lightPos = {5*256*intsin(time/4),4*65536,-4*65536+3*256*intcos(time/16)};

	// Shoot a ray through each pixel
	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
		{
			int t = 0;
			rayOrg[0] = 0;
			rayOrg[1] = 0;
			rayOrg[2] = -65536*6;

			Veci3 rayDir;
			rayDir[0] = 2*((x << 16)/XSIZE - 32768); 
			rayDir[1] = 2*((y << 16)/YSIZE*XSIZE/YSIZE/2 - 24000);
			rayDir[2] = 65536;

			#define RADIUS 384
			t = sphereIntersection(spherePos,rayOrg,rayDir,RADIUS);

			if (t < 0) 
			{
				t = planeIntersection(rayOrg,rayDir);
				if (t < 0)
				{
					t = -t;
					rayOrg[0] += (t >> 8)*(rayDir[0]) >> 4;
					rayOrg[1] += (t >> 8)*(rayDir[1]) >> 4;
					rayOrg[2] += (t >> 8)*(rayDir[2]) >> 4;

					int rx = rayOrg[0], rz = rayOrg[2];
					int dist = intsqrt((rx >> 12)*(rx >> 12)+ (rz >> 12)*(rz >> 12))/4-48;
					if (dist > 128) dist = 128;
					if (dist < 0) dist = 0;

					u8 pix = 128 - dist;

					// Calc light vector
					Veci3 L;
					L[0] = (lightPos[0] - rayOrg[0]) >> 5;
					L[1] = (lightPos[1] - rayOrg[1]) >> 5;
					L[2] = (lightPos[2] - rayOrg[2]) >> 5;

					// Shoot ray 
					t = sphereIntersection(spherePos,rayOrg,L,RADIUS);

					int lL = intsqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
					L[0] = (L[0] << 5) / lL; 
					L[1] = (L[1] << 5) / lL; 
					L[2] = (L[2] << 5) / lL; 

					if ((t < 0) || (t > lL << 4))
						setPixel(x,y,pix);
					// Point is in shade
					else
						setPixel(x,y,0);
				}
				else
				{
					setPixel(x,y,0);
				}
			} else
				{
					rayOrg[0] = rayOrg[0]+ (t*rayDir[0] >> 12);
					rayOrg[1] = rayOrg[1]+ (t*rayDir[1] >> 12);
					rayOrg[2] = rayOrg[2]+ (t*rayDir[2] >> 12);

					Veci3 N;
					N[0] = (rayOrg[0] - spherePos[0]) >> 8;
					N[1] = (rayOrg[1] - spherePos[1]) >> 8;
					N[2] = (rayOrg[2] - spherePos[2]) >> 8;

					Veci3 L;
					L[0] = (lightPos[0] - rayOrg[0]) >> 8;
					L[1] = (lightPos[1] - rayOrg[1]) >> 8;
					L[2] = (lightPos[2] - rayOrg[2]) >> 8;

					// Calculate the diffuse component of the Phong shader
					int lN = intsqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
					N[0] = (N[0] << 8) / lN; 
					N[1] = (N[1] << 8) / lN; 
					N[2] = (N[2] << 8) / lN; 

					int lL = intsqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
					L[0] = (L[0] << 8) / lL; 
					L[1] = (L[1] << 8) / lL; 
					L[2] = (L[2] << 8) / lL; 

					int angle =  (L[0]*N[0] + L[1]*N[1] + L[2]*N[2]) >> 9;

					if (angle < 0) angle = 0; 
					if (angle > 128) angle = 128;
					setPixel(x,y,angle);
				}
		}

#undef RADIUS
}

// Dithering to quantize the pixel buffer 
// to 2 color palette
// Using a sierra lite filter, see
// http://www.efg2.com/Lab/Library/ImageProcessing/DHALF.TXT
void dither()
{
	int x,y,quantError;
	for (y = 0; y < YSIZE-2; y++)
	{
		for (x = 0; x < XSIZE-2; x++)
		{
			int oldPixel = pPixelBuf[y*XSIZE+x];
			int newPixel = oldPixel/128*128;

			quantError = oldPixel - newPixel;

			pPixelBuf[y*XSIZE+x+1    ] += quantError >> 1;
			pPixelBuf[(y+1)*XSIZE+x-1] += quantError >> 2;
			pPixelBuf[(y+1)*XSIZE+x  ] += quantError >> 2;

			pPixelBuf[y*XSIZE+x] = newPixel/128*255;
		}
		int oldPixel = pPixelBuf[y*XSIZE+XSIZE-1];
		int newPixel = oldPixel/128*128;
		quantError = oldPixel - newPixel;

		pPixelBuf[(y+1)*XSIZE+XSIZE-2] += quantError >> 1;
		pPixelBuf[(y+1)*XSIZE+XSIZE-1] += quantError >> 1;
		pPixelBuf[y*XSIZE+XSIZE-1] = newPixel/128*255;
	}

	for (x = 0; x < XSIZE-2; x++)
	{
		int oldPixel = pPixelBuf[(YSIZE-1)*XSIZE+x];
		int newPixel = oldPixel/128*128;
		quantError = oldPixel - newPixel;

		pPixelBuf[(YSIZE-1)*XSIZE+x+1    ] += quantError;
		pPixelBuf[(YSIZE-1)*XSIZE+x] =  newPixel/128*255;
	}
}


void display( )
{
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);

	drawScene1();
 	//drawScene2();
	dither();

	int x,y;
	for (y = 0; y < YSIZE-1; y++)
		for (x = 0; x < XSIZE-1; x++) 
			setDispPixel(x,y,pPixelBuf[y*XSIZE+x]);

	glDrawPixels(XSIZE*SCALE-1, YSIZE*SCALE-1, GL_RGB, GL_UNSIGNED_BYTE, pDisplayBuf);

	glutSwapBuffers();
}


// Simple timer function
void update(int value)
{
	glutPostRedisplay();

	int dX = 16*intsin(time/16), dY = 256*16;
	int l = intsqrt(dX*dX + dY*dY) >> 12; 
	dX /= l; dY /= l;

	posX += dX;
	posY -= dY;
	time += 1000/24;

	glutTimerFunc(1000/24, update, 0);
}

int main(int ac, char* av[])
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
	free(pPixelBuf);

	return EXIT_SUCCESS;
}


