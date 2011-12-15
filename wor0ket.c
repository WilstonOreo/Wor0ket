#include <stdlib.h>
#include <stdio.h>
#include <GL/freeglut.h>

#define XSIZE 256
#define YSIZE 192
#define SCALE 1


#define MANT 65536.0

typedef unsigned char u8;

static u8* pPixelBuf;
static u8* pDisplayBuf;

float posX = 65536.0;
float posY = 65536.0;

typedef float Vector[3];
typedef int Veci3[3];

int lookUp_sinCos[1024];

float time = 0;

void buildLookUps()
{

}


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

static unsigned intsqrt(unsigned long x) {
	    static const unsigned char sqq_table[] = {
	       0,  16,  22,  27,  32,  35,  39,  42,  45,  48,  50,  53,  55,  57,
	      59,  61,  64,  65,  67,  69,  71,  73,  75,  76,  78,  80,  81,  83,
	      84,  86,  87,  89,  90,  91,  93,  94,  96,  97,  98,  99, 101, 102,
	     103, 104, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
	     119, 120, 121, 122, 123, 124, 125, 126, 128, 128, 129, 130, 131, 132,
	     133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 144, 145,
	     146, 147, 148, 149, 150, 150, 151, 152, 153, 154, 155, 155, 156, 157,
	     158, 159, 160, 160, 161, 162, 163, 163, 164, 165, 166, 167, 167, 168,
	     169, 170, 170, 171, 172, 173, 173, 174, 175, 176, 176, 177, 178, 178,
	     179, 180, 181, 181, 182, 183, 183, 184, 185, 185, 186, 187, 187, 188,
	     189, 189, 190, 191, 192, 192, 193, 193, 194, 195, 195, 196, 197, 197,
	     198, 199, 199, 200, 201, 201, 202, 203, 203, 204, 204, 205, 206, 206,
	     207, 208, 208, 209, 209, 210, 211, 211, 212, 212, 213, 214, 214, 215,
	     215, 216, 217, 217, 218, 218, 219, 219, 220, 221, 221, 222, 222, 223,
	     224, 224, 225, 225, 226, 226, 227, 227, 228, 229, 229, 230, 230, 231,
	     231, 232, 232, 233, 234, 234, 235, 235, 236, 236, 237, 237, 238, 238,
	     239, 240, 240, 241, 241, 242, 242, 243, 243, 244, 244, 245, 245, 246,
	     246, 247, 247, 248, 248, 249, 249, 250, 250, 251, 251, 252, 252, 253,
	     253, 254, 254, 255
	    };

	    unsigned long xn;

	    if (x >= 0x10000)
	        if (x >= 0x1000000)
	            if (x >= 0x10000000)
	                if (x >= 0x40000000) {
	                    if (x >= 65535UL*65535UL)
	                        return 65535;
	                    xn = sqq_table[x>>24] << 8;
	                } else
	                    xn = sqq_table[x>>22] << 7;
	            else
	                if (x >= 0x4000000)
	                    xn = sqq_table[x>>20] << 6;
	                else
	                    xn = sqq_table[x>>18] << 5;
	        else {
	            if (x >= 0x100000)
	                if (x >= 0x400000)
	                    xn = sqq_table[x>>16] << 4;
	                else
	                    xn = sqq_table[x>>14] << 3;
	            else
	                if (x >= 0x40000)
	                    xn = sqq_table[x>>12] << 2;
	                else
	                    xn = sqq_table[x>>10] << 1;

	            goto nr1;
	        }
	    else
	        if (x >= 0x100) {
	            if (x >= 0x1000)
	                if (x >= 0x4000)
	                    xn = (sqq_table[x>>8] >> 0) + 1;
	                else
	                    xn = (sqq_table[x>>6] >> 1) + 1;
	            else
	                if (x >= 0x400)
	                    xn = (sqq_table[x>>4] >> 2) + 1;
	                else
	                    xn = (sqq_table[x>>2] >> 3) + 1;
	
	            goto adj;
	        } else
	            return sqq_table[x] >> 4;

	/* Run two iterations of the standard convergence formula */
	    xn = (xn + 1 + x / xn) / 2;
	nr1:
	    xn = (xn + 1 + x / xn) / 2;
	adj:

	    if (xn * xn > x) /* Correct rounding if necessary */
	       xn--;
	    return xn;
	}


int sphereIntersection(Veci3 spherePos, Veci3 rayOrg, Veci3 rayDir)
{
	int t = -1000; 
	int ocVec[3];
	
	#define S 8
	#define INTRADIUS (65536)*16/10/256 

	Veci3 rayDiri;
	int i;
	for (i = 0; i < 3; i++)
	{
		ocVec[i] = (rayOrg[i] - spherePos[i]) / 256;
		rayDiri[i] = rayDir[i] / 256;
	}

	int A = rayDiri[0]*rayDiri[0] + rayDiri[1]*rayDiri[1] + rayDiri[2]*rayDiri[2];
	int B = 2*(ocVec[0]*rayDiri[0]+ocVec[1]*rayDiri[1]+ocVec[2]*rayDiri[2]);
	int C = ocVec[0]*ocVec[0]+ocVec[1]*ocVec[1]+ocVec[2]*ocVec[2]  - (INTRADIUS*INTRADIUS);

	A = A / 256; B = B / 256; C = C / 256;
	int disc = B * B - 4 * A * C;


#define BORDER 16384
//	if (disc >= -BORDER && disc <= BORDER) return 0;
	if (disc < 0)  return t;

	int distSqrt = intsqrt(disc);

	int q;
	if (B < 0)
		q = (-B - distSqrt) / 2;
	else
		q = (-B + distSqrt) / 2;

	int t0 = (!A) ? 1 << 30 : (q * 256) / A;
    int t1 = (!q) ? 1 << 30 : (C * 256) / q;

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
	if (abs(rayDir[1]) < 4) return 1 << 30; 
	else 
	 return -(((65536*4-rayOrg[1]) << 13)/rayDir[1]) >> 1;
}

u8 checkerBoard(int t, int x, int z)
{
	int dir = (t > 0) ? 1 : -1;

	int pX = (int)(posX * 256.0),pY = (int)(posY* 256.0);
	int xC = ((x*48 >> 8) + dir*pX*256) >> 16;
	int zC = ((z*48 >> 8) - dir*pY*256) >> 16;

	int checker = ((u8)xC & 1) ^ ((u8)zC & 1);

//	printf("%d %d",x,z);


	int dist = intsqrt((x >> 12)*(x >> 12)+ (z >> 12)*(z >> 12))/4-32;
	if (dist > 128) dist = 128;
	if (dist < 0) dist = 0;

	return (128-dist)*checker;
}


void draw()
{
	int x,y;
	Veci3 rayOrg;

	Veci3 spherePos;

	spherePos[1] = MANT*0.0f;
	spherePos[2] = MANT*2*cos(time/1600)+1;
	spherePos[0] = MANT*1.2*sin(time/400);
/*
	spherePos[0] = 0;
	spherePos[1] = 0;
	spherePos[2] = 0;
*/
	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++) setPixel(x,y,0);


	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
		{
			int t = 0;
			rayOrg[0] = 0;
			rayOrg[1] = 0;
			rayOrg[2] = -65536*9/2;

			Veci3 rayDir;
			rayDir[0] = 2*(x*65536/XSIZE - 32768); 
			rayDir[1] = 2*(y*65536/YSIZE*XSIZE/YSIZE/2 - 24000   );
			rayDir[2] = 65536;

			t = sphereIntersection(spherePos,rayOrg,rayDir);

			if (t < 0) 
			{
				t = planeIntersection(rayOrg,rayDir);
				rayOrg[0] += (t >> 8)*(rayDir[0]) >> 4;
				rayOrg[2] += (t >> 8)*(rayDir[2]) >> 4;
				setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
			} else
			if (t == 0)
			{
				setPixel(x,y,128);
			}
			else
			{
				rayOrg[0] += t*rayDir[0] / 256;
				rayOrg[1] += t*rayDir[1] / 256;
				rayOrg[2] += t*rayDir[2] / 256;

				Veci3 ocVec;
				ocVec[0] = (rayOrg[0] - spherePos[0]) / 256;
				ocVec[1] = (rayOrg[1] - spherePos[1]) / 256;
				ocVec[2] = (rayOrg[2] - spherePos[2]) / 256;

				int reflProd = -2*(ocVec[0]*rayDir[0] + ocVec[1]*rayDir[1] + ocVec[2]*rayDir[2])/256;
				rayDir[0] += reflProd * ocVec[0] /256 ;
				rayDir[1] += reflProd * ocVec[1] /256;
				rayDir[2] += reflProd * ocVec[2] /256;

				t = planeIntersection(rayOrg,rayDir);
				rayOrg[0] += (t >> 8)*(rayDir[0]) >> 4;
				rayOrg[2] += (t >> 8)*(rayDir[2]) >> 4;
				setPixel(x,y,checkerBoard(t,rayOrg[0],rayOrg[2])); 
			}
		}
}

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

	draw();
	dither();

	int x,y;
	for (y = 0; y < YSIZE-1; y++)
		for (x = 0; x < XSIZE-1; x++) 
			setDispPixel(x,y,pPixelBuf[y*XSIZE+x]);

	
	glDrawPixels(XSIZE*SCALE-1, YSIZE*SCALE-1, GL_RGB, GL_UNSIGNED_BYTE, pDisplayBuf);

	glutSwapBuffers();
}

void update(int value)
{
	glutPostRedisplay();

	float dX = 0.07*sin(time/2000), dY = 0.05;
	float l = sqrt(dX*dX + dY*dY); 
	dX /= l*12; dY /= l*12;

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

	return EXIT_SUCCESS;
}


