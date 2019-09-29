#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>

typedef struct {
	double real;
	double imag;
}COMPLEX;



COMPLEX com_conj(COMPLEX num) {
  COMPLEX ans;
  ans.real = num.real;
  ans.imag = -num.imag;
  return ans;
}

COMPLEX com_mul(COMPLEX num1, COMPLEX num2) {
	COMPLEX ans;
	ans.real = (num1.real*num2.real) - (num1.imag*num2.imag);
	ans.imag = (num1.real*num2.imag) + (num2.real*num1.imag);
	return ans;
}

COMPLEX com_div(COMPLEX num1, COMPLEX num2) {
  COMPLEX ans;
  ans.real = ((num1.real*num2.real)+(num1.imag*num2.imag)) / (pow(num2.real,2)+pow(num2.imag,2));
  ans.imag = ((num1.imag*num2.real)-(num1.real*num2.imag)) / (pow(num2.real,2)+pow(num2.imag,2));
  return ans;
}

COMPLEX com_addi(COMPLEX num1, COMPLEX num2) {
	COMPLEX ans = {0.0, 0.0};
	ans.real = num1.real + num2.real;
	ans.imag = num1.imag + num2.imag;
	return ans;
}

COMPLEX com_sub(COMPLEX z1,COMPLEX z2)
{
	COMPLEX ans;
	ans.real = z1.real - z2.real;
	ans.imag = z1.imag - z2.imag;
	return ans;
}

/*

COMPLEX com_div1(COMPLEX num1, COMPLEX num2) {
  COMPLEX ans;
  //ans.real = ((num1.real*num2.real)+(num1.imag*num2.imag)) / (pow(num2.real,2)+pow(num2.imag,2));
 //  return ans;
}

*/

double com_mag(COMPLEX num) {

  return sqrt(pow(num.real, 2) + pow(num.imag, 2));

}


void com_mat_mul(int ra, int ca, COMPLEX a[][ca], int rb, int cb, COMPLEX b[][cb], COMPLEX mul[ra][cb]) {
int c,d,k;
COMPLEX sum={0.0,0.0};
for (c = 0; c < ra; c++) {
	for (d = 0; d < cb; d++) {
		for (k = 0; k < rb; k++) {
			sum = com_addi(sum, com_mul(a[c][k], b[k][d]));
		}
		mul[c][d] = sum;
		sum.real = 0;
		sum.imag = 0;
	}
}
}


void com_transpose(int n, COMPLEX buff[n][n]){
int i,j;
for(i = 0; i < n; i++){
	for(j = i+1; j < n; j++){
		COMPLEX temp = buff[i][j];
		buff[i][j] = buff[j][i];
		buff[j][i] = temp;
		}
	}
}






void _fft(COMPLEX buff[], COMPLEX out[], int N, int step){
  if (step < N){
    _fft(out, buff, N, step * 2);
    _fft(out + step, buff + step, N, step * 2);

    for(int i = 0; i < N; i += 2 * step) {
      COMPLEX t = {0.0,0.0};
      t.real = cos(-i * M_PI / N);
      t.imag = sin(-i * M_PI / N);
      COMPLEX temp = com_mul(t, out[i+step]);
      buff[i/2] = com_addi(out[i], temp);
      buff[(i+N)/2] = com_sub(out[i], temp);
    }
  }
}

void fft(COMPLEX buf[], int n) {
  COMPLEX out[n];
  for (int i = 0; i < n; i++) out[i] = buf[i];
 
  _fft(buf, out, n, 1);
}

void fft_mat(int N, COMPLEX arr[N][N], COMPLEX fft_ans[N][N]) {
  int i,j;
  COMPLEX buff[N];
  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      buff[j].real = (arr[i][j].real);
      buff[j].imag = (arr[i][j].imag);
    }
    fft(buff, N);
    for(j=0; j<N; j++) {
      fft_ans[i][j].real = buff[j].real;
      fft_ans[i][j].imag = buff[j].imag;
    }
  }

  com_transpose(N, fft_ans);

  for(i=0; i<N; i++) {
    for(j=0; j<N; j++) {
      buff[j].real = fft_ans[i][j].real;
      buff[j].imag = fft_ans[i][j].imag;
    }
    fft(buff, N);
    for(j=0; j<N; j++) {
      fft_ans[i][j].real = buff[j].real;
      fft_ans[i][j].imag = buff[j].imag;
    }
  }

  com_transpose(N, fft_ans);
}

void ifft_mat(int N, COMPLEX fft[N][N], double arr[N][N]) {
  int i,j;
  COMPLEX orig[N][N];

  for(i=0;i<N;i++) {
    for(j=0; j<N; j++){
      fft[i][j] = com_conj(fft[i][j]);
    }
  }
 
  fft_mat(N, fft, orig);

  for(i=0;i<N;i++) {
    for(j=0; j<N; j++){
      arr[i][j] = com_conj(orig[i][j]).real;
      arr[i][j] /= N*N;
    }
  }
}






int main(int argc, char *argv[])
{
int N;
char in_image_fname[100], out_image_fname[100];

strcpy(in_image_fname,argv[1]);
N = atoi(argv[2]);

unsigned char (*image_buff)[N] = malloc(sizeof(unsigned char[N][N]));

int fp1,fp2;

fp1 = open(in_image_fname, O_RDONLY);

sprintf(out_image_fname,"OUTPUT.raw",in_image_fname);

fp2 = creat(out_image_fname,0667);

read(fp1, &image_buff[0][0], N*N*sizeof(unsigned char));

COMPLEX (*fft_buff)[N] = malloc(sizeof(COMPLEX[N][N]));
COMPLEX (*fft_ans)[N] = malloc(sizeof(COMPLEX[N][N]));;

for(int i=0; i<N; i++) {
	for(int j=0; j<N; j++) {
		fft_buff[i][j].real = image_buff[i][j];
		fft_buff[i][j].imag = 0;
	}
}

fft_mat(N, fft_buff, fft_ans);

COMPLEX (*psf)[N] = malloc(sizeof(COMPLEX[N][N]));
COMPLEX (*psf_ft)[N] = malloc(sizeof(COMPLEX[N][N]));

int blur = 5;
for(int i=0; i<N; i++) {
	for(int j=0; j<N; j++) {
	if(i<blur && j<blur){
		psf[i][j].real = 1/pow(blur,2);
		psf[i][j].imag = 0;
	}else {
		psf[i][j].real = 0;
		psf[i][j].imag = 0;
	}
	}
}

fft_mat(N, psf, psf_ft);

double K = 0.02;
COMPLEX (*rec_ft)[N] = malloc(sizeof(COMPLEX[N][N]));

for(int i=0; i<N; i++) {
	for(int j=0; j<N; j++) {
	COMPLEX num = com_conj(psf_ft[i][j]);
	COMPLEX den = {pow(com_mag(psf_ft[i][j]), 2) + K, 0};
	rec_ft[i][j] = com_mul(com_div(num, den), fft_ans[i][j]);
	}
}

double (*rec_img)[N] = malloc(sizeof(double[N][N]));

ifft_mat(N, rec_ft, rec_img);

for(int i=0; i<N; i++) {
	for(int j=0; j<N; j++) {
		image_buff[i][j] = (abs(round(rec_img[i][j])) > 255) ? 255 : (unsigned char)abs(round(rec_img[i][j]));
	}
}


write(fp2, &image_buff[0][0], N*N*sizeof(unsigned char));

close(fp1);
close(fp2);
return 0;
}
