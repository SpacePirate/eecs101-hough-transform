#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ROWS 480
#define COLS 640
#define DEGREES(radians) ((radians) * 57.29577951308232087679815481410517033240)
#define RADIANS(degrees) ((degrees) * 0.017453292519943295769236907684886127134)

void clear(unsigned char image[][COLS]);
void read_image(char *fname, unsigned char image[][COLS]);
void DIP(unsigned char image[][COLS], char *fname, unsigned char t);
void threshold(unsigned char image[][COLS], unsigned char t);
void header( int row, int col, unsigned char head[32] );
void save_image(char *fname, unsigned char image[][COLS]);
int find_edge(int x, int y, int rho, int theta);

int main()
{
    unsigned char image[ROWS][COLS];
    clear(image);

    read_image("image.raw", image);
    DIP(image, "image", 100);
    clear(image);

    return 0;
}

void clear( unsigned char image[][COLS] )
{
    int    i,j;
    for ( i = 0 ; i < ROWS ; i++ )
        for ( j = 0 ; j < COLS ; j++ ) image[i][j] = 0;
}

void read_image(char *fname, unsigned char image[][COLS])
{
    int i;
    FILE  *fp;
    
    /* Read in a raw image */
    /* Open the file */
    if (( fp = fopen( fname, "rb" )) == NULL )
    {
        fprintf( stderr, "error: couldn't open %s\n", fname );
        exit( 1 );
    }          
    
    /* Read the file */
    for ( i = 0; i < ROWS ; i++ )
        if ( fread( image[i], 1, COLS, fp ) != COLS )
        {
            fprintf( stderr, "error: couldn't read enough stuff\n" );
            exit( 1 );
        }

    /* Close the file */
    fclose( fp );
    
    printf("Image %s read successfully!\n", fname);
}

void save_image(char *fname, unsigned char image[][COLS])
{
    int i;
    unsigned char head[32];
    FILE  *fp;
    
    /* Save it into a ras image */
    /* Open the file */
    if (!( fp = fopen( fname, "wb" )))
      fprintf( stderr, "error: could not open %s\n", fname ), exit(1);
      
    /* Create a header */
    header(ROWS, COLS, head);
    
    /* Write the header */
    fwrite( head, 4, 8, fp );
     
    /* Write the image */
    for ( i = 0; i < ROWS; i++ )
        fwrite( image[i], 1, COLS, fp );
    
    /* Close the file */
    fclose( fp );

}
void DIP(unsigned char image[][COLS], char *fname, unsigned char t)
{
    int x, y, i, j, k;
    int sobelx, sobely, sgmval;
    int xmax = 0, ymax = 0, sgm_max = 0;
    int rho, theta;
    int minr = 0, maxr = 0;

    unsigned char edgex[ROWS][COLS];
    unsigned char edgey[ROWS][COLS];
    unsigned char SGM[ROWS][COLS];
    unsigned char hough[ROWS][COLS];

    static int bins[1600][180];
    char exname[100], eyname[100], sgmname[100], bname[100];

    clear(edgex);
    clear(edgey);
    clear(SGM);
    clear(hough);
    
    strcpy(exname, fname);
    strcpy(eyname, fname);
    strcpy(sgmname, fname);    
    strcpy(bname, fname);

    strcat(exname, "_dx.ras");
    strcat(eyname, "_dy.ras");
    strcat(sgmname, "_sgm.ras");
    strcat(bname, "_binary.ras");
    
    for(y = 1; y < ROWS-1; y++)
    {
        for(x = 1; x < COLS-1; x++)
        {
            /* dE/dx */
            sobelx = (-(image[y+1][x-1] + 2*image[y][x-1] + image[y-1][x-1]) + 
                        image[y+1][x+1] + 2*image[y][x+1] + image[y-1][x+1]);
            edgex[y][x] = sobelx;
            if(sobelx > xmax)
                { xmax = sobelx; }
                
            /* dE/dy */
            sobely = (-(image[y+1][x-1] + 2*image[y+1][x] + image[y+1][x+1]) + 
                        image[y-1][x-1] + 2*image[y-1][x] + image[y-1][x+1]);
            edgey[y][x] = sobely;
            if(sobely > ymax)
                { ymax = sobely; }
            
            /* SGM */
            sgmval = pow(pow(edgex[y][x], 2) + pow(edgey[y][x], 2), 0.5);
            SGM[y][x] = sgmval;
            if(sgmval > sgm_max)
                { sgm_max = sgmval;}
        }
    }

    /* Normalize */
    for(y = 0; y < ROWS; y++)
    {
        for(x = 0; x < COLS; x++)
        {
            edgex[y][x] = 255*edgex[y][x]/xmax;
            edgey[y][x] = 255*edgey[y][x]/ymax;
            SGM[y][x] = 255*SGM[y][x]/sgm_max;
        }
    }

    save_image(exname, edgex);
    save_image(eyname, edgey);
    save_image(sgmname, SGM);
    threshold(SGM, t);
    save_image(bname, SGM);

    /* Hough transform */
    for(y = 0; y < ROWS; y++)
    {
        for(x = 0; x < COLS; x++)
        {
            for(theta = 0; theta <= 180; theta++)
            {
                rho = -(x*sin(RADIANS(theta)) - y*cos(RADIANS(theta)));
                if (rho > maxr) {maxr = rho;}
                if (rho < minr) {minr = rho;}
                if (SGM[y][x] == 255)
                {
/*                     printf("For [%d, %d] and theta, rho is: %d, %d\n",
                            x, y, theta, rho); */
                    bins[rho + 800][theta]++;
                }
            }
        }
    }
    printf("Rho ranges from [%d, %d]\n", minr, maxr);

    /* Threshold bins */
    for(j = 0; j < 1600; j++)
    {
        for(i = 0; i < 180; i++)
        {
            if(bins[j][i] > 100)
            {
                k = j - 800;
                printf("[%d, %d] = %d\n", k, i, bins[j][i]);
            }
        }
    }

    /*
    From above, rho and theta are estimated to be:

    Line 1: rho[-289, -302], theta[126, 131], peak @ -295, 129
    Line 2: rho[-169, -171], theta[51], peak @ -171, 51
    Line 3: rho[312, 321], theta[13, 14], peak @ 312, 14
    */

    for(y = 0; y < ROWS; y++)
    {
        for(x = 0; x < COLS; x++)
        {
            /* Line 1 */
            if(find_edge(x, y, -295, 129))
                {hough[y][x]=255;}
            
            /* Line 2 */
            if(find_edge(x, y, -171, 51))
                {hough[y][x]=255;}
            
            /* Line 3 */
            if(find_edge(x, y, 312, 14))
                {hough[y][x]=255;}
        }
    }

    save_image("hough.ras",hough);

}
int find_edge(int x, int y, int rho, int theta)
{
    int a;
    a = ((x*sin(RADIANS(theta))) - (y*cos(RADIANS(theta))) + rho);
    
    if (a < 1 && a > -1)
        {return 1;}
    else
        {return 0;}
}

void threshold(unsigned char image[][COLS], unsigned char t)
{
    int x, y;
    int a = 0;
    int cx = 0, cy = 0;

    for(y = 0; y < ROWS; y++)
    {
        for(x = 0; x < COLS; x++)
        {
            /* Thresholding */
            if (image[y][x] >= t)
            {
                image[y][x] = 255;
            }
            /* Disregard background */
            else
            {
                image[y][x] = 0;
                cx += x;
                cy += y;
                a += 1;
            }
        }
    }
    
    cx = cx/a;
    cy = cy/a;
    
    /* Re-color 5x5 center area */
    for(y = cy - 2; y < cy + 3; y++)
    {
        for(x = cx - 2; x < cx + 3; x++)
        {
            image[y][x] = 128;
        }
    }

    printf("The area of the image is: %d pixels\n", a);
    printf("The center of the image is (x, y): (%d, %d)\n", cx, cy);
}

void header( int row, int col, unsigned char head[32] )
{
    int *p = (int *)head;
    char *ch;
    int num = row * col;

    /* Choose little-endian or big-endian header depending on the machine. Don't modify this */
    /* Little-endian for PC */
    
    *p = 0x956aa659;
    *(p + 3) = 0x08000000;
    *(p + 5) = 0x01000000;
    *(p + 6) = 0x0;
    *(p + 7) = 0xf8000000;

    ch = (char*)&col;
    head[7] = *ch;
    ch ++; 
    head[6] = *ch;
    ch ++;
    head[5] = *ch;
    ch ++;
    head[4] = *ch;

    ch = (char*)&row;
    head[11] = *ch;
    ch ++; 
    head[10] = *ch;
    ch ++;
    head[9] = *ch;
    ch ++;
    head[8] = *ch;
    
    ch = (char*)&num;
    head[19] = *ch;
    ch ++; 
    head[18] = *ch;
    ch ++;
    head[17] = *ch;
    ch ++;
    head[16] = *ch;
    
/*
    // Big-endian for unix
    *p = 0x59a66a95;
    *(p + 1) = col;
    *(p + 2) = row;
    *(p + 3) = 0x8;
    *(p + 4) = num;
    *(p + 5) = 0x1;
    *(p + 6) = 0x0;
    *(p + 7) = 0xf8;
*/
}