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

void draw()
{
	#define RADIUS 1.5
	int x,y;
	float ray_org[3],a,b,c,t;

	float spherePos[3];
	spherePos[0] = 1.2*sin(time/400);
	spherePos[1] = 0.0f;
	spherePos[2] = 2*cos(time/1200);

	for (y = 0; y < YSIZE; y++)
		for (x = 0; x < XSIZE; x++)
		{
			ray_org[0] = 0;
			ray_org[1] = 0;
			ray_org[2] = -3.5;

			float ray_dir[3];
			ray_dir[0] = 1.5*((float)x/(float)(XSIZE)-0.5f);
			ray_dir[1] = 1.5*((float)y/(float)(YSIZE)*0.5f*(float)XSIZE/(float)YSIZE-0.33f);
			ray_dir[2] = -1;

			float ocVec[3];
			ocVec[0] = ray_org[0] - spherePos[0];
			ocVec[1] = ray_org[1] - spherePos[1];
			ocVec[2] = ray_org[2] - spherePos[2];


			a = ray_dir[0]*ray_dir[0] + ray_dir[1]*ray_dir[1] + 1;
			b = 2*(ocVec[0]*ray_dir[0]+ocVec[1]*ray_dir[1]+ocVec[2]*ray_dir[2]);

			c = ocVec[0]*ocVec[0]+ocVec[1]*ocVec[1]+ocVec[2]*ocVec[2]  - (RADIUS*RADIUS);

			t = 0;
    		float disc = b * b - 4 * a * c;
    
#define BORDER 0.25f
			if (disc >= -BORDER && disc <= BORDER)
			{
				setPixel(x,y,128);
				continue;
			}

		    if (disc < BORDER) { //setPixel(x,y,0); 

				t = -(4.0-ray_org[1])/(ray_dir[1]);

				ray_org[0] = ray_org[0] + t*ray_dir[0];
				ray_org[1] = ray_org[1] + t*ray_dir[1];
				ray_org[2] = ray_org[2] + t*ray_dir[2];

				float dir = (t > 0) ? 1 : -1;

				float xC = (ray_org[0]*0.15) + dir*posX;
				float zC = (ray_org[2]*0.15) - dir*posY;
				if (ray_org[0] < 0) ray_org[0] -= 1;

				int   checker = ((int)xC % 2 + (int)zC%2) % 2;

				int dist = 128 - (int)sqrt(ray_org[0]*ray_org[0]/6 + ray_org[2]*ray_org[2]/6)*12;
				if (dist < 0) dist = 0; 

				setPixel(x,y,dist*(int)checker);
				continue; 
			}

			float distSqrt = sqrtf(disc);
    		float q;
    		if (b < 0)
        		q = (-b - distSqrt)/2.0;
    		else
        		q = (-b + distSqrt)/2.0;

    		float t0 = q / a;
   			//float t1 = c / q;


    		//if (t1 < 0) continue;
				
			if (t0 < 0)
    		{
        		t = t0;
				setPixel(x,y,128); 
    		}
    		else
    		{
        		t = t0;
				setPixel(x,y,64); 
    		}

			ray_org[0] = ray_org[0] + t*ray_dir[0];
			ray_org[1] = ray_org[1] + t*ray_dir[1];
			ray_org[2] = ray_org[2]  + t*ray_dir[2];
			ocVec[0] = ray_org[0] - spherePos[0];
			ocVec[1] = ray_org[1] - spherePos[1];
			ocVec[2] = ray_org[2] - spherePos[2];


			float reflProd = 2*(ocVec[0]*ray_dir[0] + ocVec[1]*ray_dir[1] + ocVec[2]*ray_dir[2]);

			ray_dir[0] = ray_dir[0] - reflProd * ocVec[0] ;
			ray_dir[1] = ray_dir[1] - reflProd * ocVec[1] ;
			ray_dir[2] = ray_dir[2] - reflProd * ocVec[2] ;

			t = -(4.0-ray_org[1])/(ray_dir[1]);
	//		printf("%f \t%f \t%f \t --- %f \t%f \t%f \t --- %f\n",ray_org[0],ray_org[1],ray_org[2],ray_dir[0],ray_dir[1],ray_dir[2],t);

			//if (t > 0)
			{
				ray_org[0] = ray_org[0] + t*ray_dir[0];
				ray_org[1] = ray_org[1] + t*ray_dir[1];
				ray_org[2] = ray_org[2] + t*ray_dir[2];

				float dir = (t > 0) ? 1 : -1;

				float xC = (ray_org[0]*0.15) + dir*posX;
				float zC = (ray_org[2]*0.15) - dir*posY;
				if (ray_org[0] < 0) ray_org[0] -= 1;

				int   checker = ((unsigned int)xC % 2 + (unsigned int)zC % 2) % 2;
				

				int dist = 128 - (int)(ray_org[0]*ray_org[0] + ray_org[2]*ray_org[2])/3;
				if (dist < 0) dist = 0; 
				
				setPixel(x,y,(128*checker)*dist/128);
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
			//printf("%d ",newPixel/64);

			int quantError = oldPixel - newPixel;
			
			/* // Floyd-Steinberg
		    pPixelBuf[y*XSIZE+x+1    ] += 7 * quantError / 16;
     		pPixelBuf[(y+1)*XSIZE+x-1] += 3 * quantError / 16;
      	    pPixelBuf[(y+1)*XSIZE+x  ] += 5 * quantError / 16;
      		pPixelBuf[(y+1)*XSIZE+x+1] += 1 * quantError / 16;*/
		
			pPixelBuf[y*XSIZE+x+1    ] += quantError / 4;
			pPixelBuf[y*XSIZE+x+2    ] += quantError / 8;
			pPixelBuf[(y+1)*XSIZE+x-1] += quantError / 8;
			pPixelBuf[(y+1)*XSIZE+x  ] += quantError / 4;
			pPixelBuf[(y+1)*XSIZE+x+1] += quantError / 8;
			pPixelBuf[(y+2)*XSIZE+x  ] += quantError / 8;

		

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

	posX += 0.07*sin(time/2000);
	posY += 0.05;
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


