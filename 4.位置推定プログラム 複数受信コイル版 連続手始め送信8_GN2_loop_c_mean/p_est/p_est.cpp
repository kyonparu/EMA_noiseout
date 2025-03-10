#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <float.h> //FLT_MAXを使うため
#include <limits.h>
#include <string.h>
#include <direct.h>

#define NR_END 1
#define FREE_ARG char*

#define SN 6 //計算したい行列の長さよりも1大きい数をいれてね！3x3なら4
#define CL 60.0    //コイル長[m]
#define L  337.5    //座標中心からコイルの中心までの長さ[m]
#define M_PI 3.141592653589793
#define T_sample 50000*3 /*1チャンネル分のサンプル数*/
#define K 8 /*周波数分解点数*/
#define N_of_channels 12 /*チャンネル数*/
#define sect 400
#define FS 50000
#define PATERN 17
#define TEST 2187
//#define PARMIT_ERROR_LEVEL 15.0
//#define PARMIT_ERROR_LEVEL 22.0
#define PARMIT_ERROR_LEVEL 20.0

#define NUM_SENT_COIL 8 //送信コイル数
#define AXIS_NUM 5 //状態の次元数
#define SUB_AXIS_NUM 3 //位置の次元数
#define NUM_POLE 2 //極の数
#define FILENUM 1 //校正時使用ファイル数
#define DATA_FILENUM 1 //データファイル数

#define START_CH_NO 1 //最初の受信コイルNo.
#define END_CH_NO 12 //最後の受信コイルNo.
#define NUM_RES_COIL (END_CH_NO - START_CH_NO + 1) //受信コイル数
#define REF_CH_NO 2

#define SCOIL_DIR "2025022703"

#define DATA_DATE "20250304"
#define START_DATA_NO 10801
#define END_DATA_NO 10802
//#define DATA_DIR "2015110911mov2"
//#define DATA_FILE_NAME_BEFOR_FILE_NO "AD2015110911_out"

int min_index(double *S,int i1,int i2){
int i,j;
double C;
j=i1;
C=S[i1];
for(i=i1+1;i<i2;i++){
	if(S[i] < C){
		j=i;
		C=S[i];
		}
	}

return(j);
}
int max_index(double *S,int i1,int i2){
int i,j;
double C;
j=i1;
C=S[i1];
for(i=i1+1;i<i2;i++){
	if(S[i] > C){
		j=i;
		C=S[i];
		}
	}

return(j);
}
double **new_matrix(int m, int n)
{
  int i,j;
  double **a;
  a = (double **)malloc(sizeof(double *) * m);
  for (i = 0; i < m; i++)
	  a[i] = (double *)malloc(sizeof(double) * n);

for(i=0;i<m;i++){
	for(j=0;j<n;j++){
		a[i][j]=0.0;
		}
	}
  return a;
}
void free_matrix(double **a,int m)
{
	int i;

	for (i = 0; i < m; i++){
		free(a[i]);
//		printf("%dunnnko\n",i);
		a[i]=NULL;
	}

	free(a);
	a=NULL;
}
double **Data_import(char Datafilename[],char ScoilDirname[],char Scoilfilename[],char GPfilename[] ,char Spfilename[], int *Lines,double **Scoil,double *Qm,double *r, double **dd_read){
//double Scoil[][]だとダメだった。二次元配列だと、どっちか要素が要るらしい
int i,j, index;
char *inbuf;
double buf[NUM_SENT_COIL];
double temp_Qm[NUM_SENT_COIL];
double temp_Scoil[SUB_AXIS_NUM * NUM_POLE][NUM_SENT_COIL];
FILE *fp;
double **data_buf;
char flname[128];



sprintf(flname,"%s.data",Datafilename);
//ファイルオープン
if((fp = fopen(flname, "r")) == NULL ) {
	printf("ファイルオープンエラー @data file 1\n");
	exit(EXIT_FAILURE);
}
//////////////

//memory allocation
	inbuf = (char *)malloc(128);
//////////////
	i=0;
	while(fgets(inbuf,128,fp)!=NULL) i++;//これは行数を数えるだけ
	fclose(fp);
	printf("file ended %d Lines\n",i);
	*Lines=i;
// 	//// 確認用------------------------------
// 	printf("%d:Lines\n",*Lines);
// 	printf("%s:inbuf\n",inbuf);	
// 	////------------------------------------
	
	//// memory allocation
	data_buf=new_matrix(*Lines,NUM_SENT_COIL);
	////------------------

	//ファイルオープン
	if((fp = fopen(flname, "r")) == NULL ) {
		printf("ファイルオープンエラー @data file 2\n");
		exit(EXIT_FAILURE);
	}
	//////////////	
	i=0;
// 	while(fgets(inbuf,128,fp)!=NULL){
// 		printf("unnko %d\n",i);
// 		sscanf(inbuf, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf ",&data_buf[i][0],&data_buf[i][1],&data_buf[i][2],&data_buf[i][3],&data_buf[i][4],&data_buf[i][5]);
// 		i++;
// 	}
	
	while(fgets(inbuf,128,fp)!=NULL){
//		printf("unnko %d\n",i);
		sscanf(inbuf, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf ",&buf[0],&buf[1],&buf[2],&buf[3],&buf[4],&buf[5], &buf[6], &buf[7]);

		for(j=0;j<NUM_SENT_COIL;j++){
			data_buf[i][j]=buf[j];	//こういう風に一回バッファに入れてからじゃないとdata_bufには直接入れられなかった。
		}
		i++;
	}

	fclose(fp);

// 	////// 確認用--------------------
// 	for(i=0;i<*Lines;i++){
// 		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",data_buf[i][0],data_buf[i][1],data_buf[i][2],data_buf[i][3],data_buf[i][4],data_buf[i][5]);
// 	}
// 	/////----------------------------
	
	
	/////////////   ここからScoilのデータ読み込み
	sprintf(flname,"%s_Scoil.data",Scoilfilename);
// 	sprintf(flname,"../Calib_Cdouble/Results/%s/%s/%s_Scoil.txt\0",ScoilDirname,Scoilfilename,Scoilfilename);
	//ファイルオープン
	if((fp = fopen(flname, "r")) == NULL ) {
		printf("ファイルオープンエラー @Scoil\n");
		exit(EXIT_FAILURE);
	}
	/////------------
	i=0;
	while(fgets(inbuf,128,fp)!=NULL){
		sscanf(inbuf, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t ",&buf[0],&buf[1],&buf[2],&buf[3],&buf[4],&buf[5]);
		for(j=0;j<6;j++){
			Scoil[i][j]=buf[j];	//こういう風に一回バッファに入れてからじゃないとdata_bufには直接入れられなかった。
		}
		i++;
	}
	
	fclose(fp);

	/////////////   ここからQm(ゲインパラメータ)のデータ読み込み
	sprintf(flname,"%s_GP.data",GPfilename);
// 	sprintf(flname,"../Calib_Cdouble/Results/%s/%s/%s_GP.txt\0",ScoilDirname,Scoilfilename,GPfilename);
	//ファイルオープン
	if((fp = fopen(flname, "r")) == NULL ) {
		printf("ファイルオープンエラー @GP\n");
		exit(EXIT_FAILURE);
	}
	/////------------
	//while(fgets(inbuf,128,fp)!=NULL)
	//{
		fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&buf[0],&buf[1],&buf[2],&buf[3],&buf[4],&buf[5], &buf[6], &buf[7]);
		for(j=0;j<NUM_SENT_COIL;j++)
		{
			Qm[j]=buf[j];	//こういう風に一回バッファに入れてからじゃないとdata_bufには直接入れられなかった。
		}
		fscanf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",&buf[0],&buf[1],&buf[2],&buf[3],&buf[4],&buf[5], &buf[6], &buf[7]);
		for(j=0;j<NUM_SENT_COIL;j++)
		{
			r[j]=buf[j];	//こういう風に一回バッファに入れてからじゃないとdata_bufには直接入れられなかった。
		}
		printf( "G:%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Qm[0],Qm[1],Qm[2],Qm[3],Qm[4],Qm[5], Qm[6], Qm[7]);
		printf( "R:%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r[0],r[1],r[2],r[3],r[4],r[5], r[6], r[7]);
	//}

		fclose(fp);
	
// 	////// 確認用--------------------
// 	puts("Scoil_data");
// 	for(i=0;i<6;i++){
// 		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Scoil[i][0],Scoil[i][1],Scoil[i][2],Scoil[i][3],Scoil[i][4],Qm[i]);
// 	}
// 	/////----------------------------
// 	
		
		/////////////   ここからSp(受信コイル位置ずれ量)のデータ読み込み
	sprintf(flname,"%s_Sp.data", Spfilename);
// 	sprintf(flname,"../Calib_Cdouble/Results/%s/%s/%s_GP.txt\0",ScoilDirname,Scoilfilename,GPfilename);
	//ファイルオープン
	if((fp = fopen(flname, "r")) == NULL ) {
		printf("ファイルオープンエラー @GP\n");
		exit(EXIT_FAILURE);
	}
	/////------------
	//while(fgets(inbuf,128,fp)!=NULL)
	//{
	printf("test1\n");
	for(i = 0; i < NUM_SENT_COIL; i++)
		fscanf(fp, "%lf\n", &temp_Qm[i]);

	for (i = 0; i < NUM_SENT_COIL; i++) {
		for (j = 0; j < NUM_POLE * SUB_AXIS_NUM - 1; j++)
			fscanf(fp, "%lf\t", &temp_Scoil[j][i]);
		fscanf(fp, "%lf\n", &temp_Scoil[j][i]);
	}

	for (i = 0; i < FILENUM; i++) {
		for (j = 0; j < AXIS_NUM - 1; j++)
			fscanf(fp, "%lf\t", &dd_read[i][j]);
		fscanf(fp, "%lf\n", &dd_read[i][j]);
	}

	fscanf(fp, "%d\n", &index);
	printf("dd:");

	for (i = 0; i < FILENUM; i++) {
		for (j = 0; j < AXIS_NUM - 1; j++)
			printf("%f\t",dd_read[i][j]);
		printf("%f\n", dd_read[i][j]);
	}
	//}
	
// 	////// 確認用--------------------
// 	puts("Scoil_data");
// 	for(i=0;i<6;i++){
// 		printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",Scoil[i][0],Scoil[i][1],Scoil[i][2],Scoil[i][3],Scoil[i][4],Qm[i]);
// 	}
// 	/////----------------------------
// 	
	
	printf("test2\n");

	free(inbuf);

	printf("test3\n");
	fclose(fp);
	printf("test4\n");
	return(data_buf);
}

void angle_wrap1(double *a){
while(*a <= -M_PI)
	*a = *a + 2.0*M_PI;
	
while(*a > M_PI)
	*a = *a - 2.0*M_PI;
}
void angle_wrap2(double *a){
	while(*a <= -M_PI+0.6)
	*a = *a + 2.0*M_PI;
	
	while(*a > M_PI+0.6)
		*a = *a - 2.0*M_PI;
}/////表示するときはこっちね
void Dfunc(double *X,double qm,double *Scoil_pole_coord,double CC,double CS,double Si,double *Z){
//送信コイルの極がScoil_pole_coordにあり、その磁荷がqmのときに
//座標Xにあり角度はCC,CS,Siであらわされる受信コイルに予測されるゲインを計算する
// X …… 測定点（受信コイル位置）の座標(m)
// qm …… 送信コイルの磁荷
// qm=1で受け取ればミアケくんのAfuncと一緒だよ

int i;
double Q,Rp[3],Rn[3];
double lenRp,lenRn,Yp[3],Yn[3],Y[3];
lenRp=0.0;
lenRn=0.0;
//U0 = 1.26/(10^6);      //空気の透磁率≒真空の透磁率＝1.26×10^-6

	Q=qm/(4.0*M_PI*(1.26/pow(10.0,4.0)));
//	Q=qm/(16.0*M_PI*M_PI/pow(10,10));//理論どおりではたぶん10
////----データ入力確認用-----
// printf("X:%lf,%lf,%lf\n",X[0],X[1],X[2]);
// for(i=0;i<6;i++){
// 	printf("Scoil_pole_coord[%d]:%lf\n",i,Scoil_pole_coord[i]);
// }
// printf("qm:%lf\n",qm);
// printf("Q:%lf\n",Q);
////----------------------

//S:1x6 coilNo番目の両極の座標
for(i=0;i<3;i++){
    Rp[i] = X[i] - Scoil_pole_coord[i];             //p極から点X(今の受信信号の位置)までのベクトル
    Rn[i] = X[i] - Scoil_pole_coord[i+3];             //n極から点Xまでのベクトル

	lenRp=lenRp+pow(Rp[i],2.0);
	lenRn=lenRn+pow(Rn[i],2.0);
	}
lenRp=pow(lenRp,0.5);
lenRn=pow(lenRn,0.5);
////lenRp確認用---	
// printf("lenRp:%lf\n",lenRp);
// printf("lenRn:%lf\n",lenRn);
////-------------

for(i=0;i<3;i++){
	Yp[i] =  Q*(Rp[i]/pow(lenRp,3.0));   //p極の磁荷が測定点Xにつくる磁場
	Yn[i] = -Q*(Rn[i]/pow(lenRn,3.0));  //n極の磁荷が測定点Xにつくる磁場（負電荷なのでマイナス）
    Y[i] = Yp[i] + Yn[i];
    //各送信コイルの両極の磁荷が点Xにつくる磁場のベクトル（1x3、各行がコイル番号に対応）
	}
////lenRp確認用---
// for(i=0;i<3;i++){
// 	printf("Yp[%d]:%lf\n",i,Yp[i]);
// 	printf("Yn[%d]:%lf\n",i,Yn[i]);
// 	printf("Y[%d]:%lf\n",i,Y[i]);
// 	}
////-------------

Z[0] = Y[0]*CC + Y[1]*CS - Y[2]*Si;
/////CC,CS,Si確認用-----
// printf("CC:%lf\n",CC);
// printf("CS:%lf\n",CS);
// printf("Si:%lf\n",Si);
// printf("Z:%lf\n",Z[0]);
//予測ゲイン

}
void Dfuncg_dd(double *X,double qm,double *s,double *G,int T2,int I,double *bb){
int i;
double Q,Rp[3],Rn[3];
double lenRp=0,lenRn=0,Yp[3],Yn[3],Y[3],CC2,CS2,Si2,CC,CS,Si;
lenRp=0.0;
lenRn=0.0;
//U0 = 1.26/(10^6);      //空気の透磁率≒真空の透磁率＝1.26×10^-6

CC = cos(bb[3]/180.0*M_PI)*cos(bb[4]/180.0*M_PI);
CS = cos(bb[3]/180.0*M_PI)*sin(bb[4]/180.0*M_PI);
Si = -sin(bb[3]/180.0*M_PI);


CC2=CC;
CS2=CS;
Si2=Si;


	
    Q=qm/(4.0*M_PI*(1.26/pow(10.0,4.0)));
//    Q=qm/(4.0*M_PI*(1.26/pow(10,4)));
//printf("X1:%f Y1:%f Z1:%f X2:%f Y2:%f Z2:%f\n",s[0],s[1],s[2],s[3],s[4],s[5]);
//printf("%f %f %f\n",X[T2][0],X[T2][1],X[T2][2]);
//	printf("dd%lf\t%lf\t%lf\t%lf\t%lf\n",bb[0],bb[1],bb[2],bb[3],bb[4]);
//	printf("test%d\n",DD);

	//dd[0]=0;dd[1]=0;dd[2]=0;
//S:1x6 coilNo番目の両極の座標

    Rp[0] = X[0] + bb[0] - s[0];             //p極から点X(今の受信信号の位置)までのベクトル
    Rp[1] = X[1] + bb[1] - s[1];
    Rp[2] = X[2] + bb[2] - s[2];
    Rn[0] = X[0] + bb[0] - s[3];             
    Rn[1] = X[1] + bb[1] - s[4];
    Rn[2] = X[2] + bb[2] - s[5];            //n極から点Xまでのベクトル


for(i=0;i<3;i++){
//printf("Rp[%d]:%lf ",i,Rp[i]);
//printf("Rn[%d]:%lf\n",i,Rn[i]);
    lenRp=lenRp+pow(Rp[i],2);
    lenRn=lenRn+pow(Rn[i],2);

}
//puts("ok");
    
lenRp=pow(lenRp,0.5);
lenRn=pow(lenRn,0.5);
////lenRp確認用---    
//printf("lenRp %lf\n",lenRp);
//printf("lenRn %lf\n",lenRn);
////-------------

for(i=0;i<3;i++){
    Yp[i] =  Q*(Rp[i]/pow(lenRp,3));   //p極の磁荷が測定点Xにつくる磁場
    Yn[i] = -Q*(Rn[i]/pow(lenRn,3));  //n極の磁荷が測定点Xにつくる磁場（負電荷なのでマイナス）
    Y[i] = Yp[i] + Yn[i];
    //各送信コイルの両極の磁荷が点Xにつくる磁場のベクトル（1x3、各行がコイル番号に対応）
    }
////lenRp確認用---
for(i=0;i<3;i++){
//    printf("Yp[%d]:%lf ",i,Yp[i]);
//   printf("Yn[%d]:%lf ",i,Yn[i]);
//   printf("Y[%d]:%lf\n",i,Y[i]);
    }
////-------------

G[T2] = Y[0]*CC2 + Y[1]*CS2 + Y[2]*Si2;

//if(G[c]<0.0)G[c]=-G[c];

//if(Z>0.00053 || Z< -0.00053)return(1);else return(0);

//GS[H][I][S][T2][c]=G[c];




/////CC,CS,Si確認用-----
//    printf("CC:%lf\n",CC);
//    printf("CS:%lf\n",CS);
//    printf("Si:%lf\n",Si);
//    printf("Z[%d][%d] %lf\n",c,T2,G[T2]);
//予測ゲイン



}
double *vector(long nl, long nh)
{
 double *v;

 v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
 if(!v) puts("allocation error in vector()");
 return v-nl+1;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */ 
{
	int *v;

	v = (int *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(int)));

	if (!v)
		puts("allocation failure in ivector()");

	return v - nl + NR_END;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG)(v + nl - NR_END));
}


int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
	int **m;

	/* allocate pointers to rows */
	m = (int **) malloc((size_t)((nrow + NR_END) * sizeof(int *)));

	if (!m) 
		puts("allocation failure 1 in imatrix()");

	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl] = (int *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));

	if (!m[nrl])
		puts("allocation failure 2 in imatrix()");

	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for (i = nrl + 1; i <= nrh; i++)
		m[i] = m[i - 1] + ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG)(m[nrl] + ncl - NR_END));
	free((FREE_ARG)(m + nrl - NR_END));
}

void free_vector(double *v, long nl, long nh)
{
 free( (char*) (v+nl-1) );
}
int ludcmp(double **a, int n, int *indx, double *d)
{
 int i,imax,j,k;
 double big,dum,sum,temp;
 double *vv; //各行の暗黙のスケーリングの記録

 vv=vector(1,n);//ここをどうにかしないと


 *d=1.0;

 for(i=1;i<=n;i++){//行のループ
  big=0.0;
  for(j=1;j<=n;j++)
   if((temp=fabs(a[i][j]))>big) big=temp;
  if(big ==0.0) {
	  puts("Singular matrix in routine ludcmp 特異行列です");//最大要素が０なら特異行列
//	  return(1);
  }
  vv[i]=1.0/big;//スケーリング記録
 }

 for(j=1;j<=n;j++){ //Crout法の列のループ
  for(i=1;i<j;i++){
   sum=a[i][j];
   for(k=1;k<i;k++) sum -= a[i][k]*a[k][j];
   a[i][j]=sum;
  }
  big=0.0;
  for(i=j;i<=n;i++){
   sum=a[i][j];
   for(k=1;k<j;k++)
    sum -= a[i][k]*a[k][j];
   a[i][j]=sum;
   if( (dum=vv[i]*fabs(sum))>=big){
    big=dum;
    imax=i;
   }
  }

  if(j != imax){
   for(k=1;k<=n;k++){
//	puts("unnkoX");
    dum=a[imax][k];
//	puts("unnkoX2");
    a[imax][k]=a[j][k];
    a[j][k]=dum;
   }

   *d=-(*d);
   vv[imax]=vv[j];

  }

  indx[j]=imax;
  
  if(a[j][j]==0.0) a[j][j]=1.0e-20;

  if(j!=n){
   dum=1.0/(a[j][j]);
   for(i=j+1;i<=n;i++) a[i][j] *=dum;
  }
 }


 free_vector(vv,1,n);
 return(0);
}
void lubksb(double **a, int n, int *indx, double b[])
{
 int i,ii=0,ip,j;
 double sum;

 for(i=1;i<=n;i++){
  ip=indx[i];
  sum=b[ip];
  b[ip]=b[i];
  if(ii)
   for(j=ii;j<=i-1;j++) sum -=a[i][j]*b[j];
  else if(sum) ii=i;
  b[i]=sum;
 }

 for(i=n;i>=1;i--){
  sum=b[i];
  for(j=i+1;j<=n;j++) sum -=a[i][j]*b[j];
  b[i]=sum/a[i][i];
 }
}

void Calc_Inv_main(double **Ain, double **Yout)
{
double **a, d,**y,*col;
int p,n,*indx;
int i,j;




//領域確保
col=(double *)malloc(sizeof(double)*SN);
indx=(int *)malloc(sizeof(int)*SN);
//b=(double *)malloc(sizeof(double)*SN);



a=new_matrix(SN,SN);
y=new_matrix(SN,SN);
//Yout=new_matrix(SN-1,SN-1);


///////======入力の行列がインデックス0からのに対応させるため=====
for(i=0;i<SN-1;i++){
 for(j=0;j<SN-1;j++){
	a[i+1][j+1]=Ain[i][j];
	}
 }////塚本君のは式とインデックスを一致させるためにインデックス0を使ってないから
/////============================================================

// n は扱う行列の大きさです　SNはただの#define
n=SN-1;
d=0.0;

//イニシャライズ
for(p=0;p<=n;p++)
indx[p]=0;


//ここから逆行列のためのくだり

if(ludcmp(a,n,indx,&d)==1)
	return;



for(j=1;j<=n;j++){
 for(i=1;i<=n;i++)col[i]=0.0;
  col[j]=1.0;
  lubksb(a,n,indx,col);
  for(i=1;i<=n;i++) y[i][j]=col[i];
}

///////----出力するためにインデックスを0からにした。
for(i=0;i<SN-1;i++){
 for(j=0;j<SN-1;j++){
	Yout[i][j]=y[i+1][j+1];
	}
 }
///////-------------------------------------------

free(col);
free(indx);
//free(b);

free_matrix(a,SN);
free_matrix(y,SN);


//return(Yout);
}

//double **Calc_Inv_main(double **Ain)
//{
//double **a, d,**y,*col,**Yout;
//int p,n,*indx;
//int i,j;
//
//
//
//
////領域確保
//col=(double *)malloc(sizeof(double)*SN);
//indx=(int *)malloc(sizeof(int)*SN);
////b=(double *)malloc(sizeof(double)*SN);
//
//
//
//a=new_matrix(SN,SN);
//y=new_matrix(SN,SN);
//Yout=new_matrix(SN-1,SN-1);
//
//
/////////======入力の行列がインデックス0からのに対応させるため=====
//for(i=0;i<SN-1;i++){
// for(j=0;j<SN-1;j++){
//	a[i+1][j+1]=Ain[i][j];
//	}
// }////塚本君のは式とインデックスを一致させるためにインデックス0を使ってないから
///////============================================================
//
//// n は扱う行列の大きさです　SNはただの#define
//n=SN-1;
//d=0.0;
//
////イニシャライズ
//for(p=0;p<=n;p++)
//indx[p]=0;
//
//
////ここから逆行列のためのくだり
//
//if(ludcmp(a,n,indx,&d)==1)return(NULL);
//
//
//
//for(j=1;j<=n;j++){
// for(i=1;i<=n;i++)col[i]=0.0;
//  col[j]=1.0;
//  lubksb(a,n,indx,col);
//  for(i=1;i<=n;i++) y[i][j]=col[i];
//}
//
/////////----出力するためにインデックスを0からにした。
//for(i=0;i<SN-1;i++){
// for(j=0;j<SN-1;j++){
//	Yout[i][j]=y[i+1][j+1];
//	}
// }
/////////-------------------------------------------
//
//free(col);
//free(indx);
////free(b);
//
//free_matrix(a,SN);
//free_matrix(y,SN);
//
//
//return(Yout);
//}
int Calc_dd(double **Zy,int RLPi,double **POS,double **angles,int a,double **Scoil_pole,double *Qm,double *S_return){
// 更新ベクトル Δdを計算してそれを加え、更新した受信コイル推定位置と角度を返す。
// また、そのときの予測誤差規範Sも返す
// 引数は実測ゲイン:[[0]]x6  受信コイル位置、角度、送信コイルの位置角度と磁荷6x6
int i,j,k,I;
double Q[NUM_SENT_COIL];
double estG0[NUM_SENT_COIL],dd[5];
double CC,CS,Si,SC,SS,Ci;

double Rp[3],Rn[3],lenRp,lenRn;
double roRp3x1,roRp3x2,roRp3x3,roRn3x1,roRn3x2,roRn3x3;
double **A,**B,**C,**InvC,**InvCxB;
double **Gdiff;
double M[]={5.0,1.0,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005, 0.00001, 5.0e-6, 1.0e-6, 5.0e-7, 1.0e-7, 5.0e-8, 1.0e-8};//更新ベクトルにかける重み
int Mlen =18;
// double M[]={1};
// int Mlen=1;
double **POSbuff,**anglebuff;
double S[18];




A=new_matrix(NUM_SENT_COIL,5);
B=new_matrix(5,NUM_SENT_COIL);
C=new_matrix(5,5);
InvC=new_matrix(SN-1,SN-1);
//InvCxB=new_matrix(5,5);
InvCxB = new_matrix(5, NUM_SENT_COIL);
Gdiff=new_matrix(NUM_SENT_COIL,1);
POSbuff=new_matrix(Mlen,3);
anglebuff=new_matrix(Mlen,2);

for(j=0;j<NUM_SENT_COIL;j++){
//printf("%d,Q[j]%lf \n",j,Qm[j]);
	Q[j]=Qm[j]/(4.0*M_PI*(1.26/pow(10.0,4.0)));
//	Q[j]=Q[j]/(16.0*M_PI*M_PI/pow(10,10));//理論どおりではたぶん10
}

//-----------------------ヤコビ行列算出-------by Okamura--------------------//
////---まず今の受信コイル位置での予測信号ゲインを計算する それを位置情報で偏微分する

CC = cos(angles[a][0])*cos(angles[a][1]);
CS = cos(angles[a][0])*sin(angles[a][1]);
Si = sin(angles[a][0]);  //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある
SC = sin(angles[a][0])*cos(angles[a][1]);
SS = sin(angles[a][0])*sin(angles[a][1]);
Ci = cos(angles[a][0]);


for(i=0;i<NUM_SENT_COIL;i++){ //i番目の送信コイルに対して
//printf("Q[j]%lf \n",Q[j]);	
	for(j=0;j<3;j++){
	Rp[j]=0.0;
	Rn[j]=0.0;
	Rp[j] = POS[a][j] - Scoil_pole[i][j];//i番目送信コイル正極から受信コイル位置までのベクトル
    Rn[j] = POS[a][j] - Scoil_pole[i][j+3];//
	}

lenRp=0.0;
lenRn=0.0;

	for(j=0;j<3;j++){
	lenRp=lenRp+pow(Rp[j],2.0);
	lenRn=lenRn+pow(Rn[j],2.0);   //ベクトルの長さ
	}

	lenRp=pow(lenRp,0.5);
	lenRn=pow(lenRn,0.5);


    roRp3x1 = 3.0 * lenRp * Rp[0];  //偏微分らしいけどなんでこうなるの？やまさんの論文9ページ
    roRp3x2 = 3.0 * lenRp * Rp[1];
    roRp3x3 = 3.0 * lenRp * Rp[2];
    roRn3x1 = 3.0 * lenRn * Rn[0];
    roRn3x2 = 3.0 * lenRn * Rn[1];
    roRn3x3 = 3.0 * lenRn * Rn[2];

  A[i][0] = Q[i] * (  ( (pow(lenRp,3) - Rp[0]*roRp3x1) / (pow(lenRp,6)) - (pow(lenRn,3) - Rn[0]*roRn3x1) / (pow(lenRn,6)) ) * CC           + ( ((-Rp[1])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x1) / (pow(lenRn,6)) ) * CS
                  - ( ((-Rp[2])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x1) / (pow(lenRn,6)) ) * Si  );


  A[i][1] = Q[i] * (  ( ((-Rp[0])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x2) / (pow(lenRn,6)) ) * CC
           +( (pow(lenRp,3) - (Rp[1]*roRp3x2)) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[1])*roRn3x2) / (pow(lenRn,6)) ) * CS
                 - ( ((-Rp[2])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x2) / (pow(lenRn,6)) ) * Si);
             
  A[i][2] = Q[i] * (  ( ((-Rp[0])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x3) / (pow(lenRn,6)) ) * CC
           +( ((-Rp[1])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x3) / (pow(lenRn,6)) ) * CS
                - ( (pow(lenRp,3) - (Rp[2])*roRp3x3) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[2])*roRn3x3) / (pow(lenRn,6)) ) * Si  );

  A[i][3] = Q[i] * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*SC
                 -( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*SS
                  -( (Rp[2])/(pow(lenRp,3)) - (Rn[2])/(pow(lenRn,3)) )*Ci   );
                  
  A[i][4] = Q[i] * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*CS
                 +( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*CC   );  

 //Aはヤコビアン行列（6x5、各行がコイル番号）
//printf("Q[i]%lf \n",Q[i]);

//printf("%lf %lf %lf %lf %lf\n",A[i][0],A[i][1],A[i][2],A[i][3],A[i][4]);
}
//--------------------------------------------------------------------//
//---------------------------受信信号予測値算出-------------------------//



for(i=0;i<NUM_SENT_COIL;i++){
Dfunc(&POS[a][0],Qm[i],&Scoil_pole[i][0],CC,CS,Si,&estG0[i])  ;  //受信信号予測値（NUM_SENT_COILx[[0]]、各行が送信コイル番号に対応）;
}
//--------------------------------------------------------------------//

//printf("B[i][j] \n");
////----転地行列-Bの計算---------------
for(i=0;i<5;i++){
	for(j=0;j<NUM_SENT_COIL;j++){
		B[i][j]=A[j][i];
	}
//printf("%lf %lf %lf %lf %lf %lf\n",B[i][0],B[i][1],B[i][2],B[i][3],B[i][4],B[i][5]);
}
//----------B:5x6-------------------


////=====行列掛け算  C 5x5 ========
//printf("C[i][j] \n");
for(i=0;i<5;i++){
	for(j=0;j<5;j++){
		for(k=0;k<NUM_SENT_COIL;k++){
			C[i][j]+=B[i][k]*A[k][j];
		}
	}

//printf("%lf %lf %lf %lf %lf\n",C[i][0],C[i][1],C[i][2],C[i][3],C[i][4]);
}//==============================



//////--逆行列の計算------------------------
Calc_Inv_main(C, InvC);
//InvC=Calc_Inv_main(C); // Calc_Inv_mainは塚本氏によるLU分解を用いた逆行列計算、2は吐き出し法を用いた逆行列計算
//if(InvC==NULL){
//printf("InvC=NULL!\n");
//return((int)DBL_MAX);
//}
//////------------------------------------
//puts("unnko4");
////=====行列掛け算  InvCxB 5x5 ==========================
//for(i=0;i<5;i++){
//	for(j=0;j<5;j++){
//		for(k=0;k<5;k++){
//			InvCxB[i][j]+=InvC[i][k]*B[k][j];
//		}
//	}
//}
for(i=0;i<5;i++){
	for(j=0;j<NUM_SENT_COIL;j++){
		for(k=0;k<5;k++){
			InvCxB[i][j]+=InvC[i][k]*B[k][j];
		}
	}
}
//=====================================================

for(j=0;j<NUM_SENT_COIL;j++){
	Gdiff[j][0]=Zy[RLPi][j]-estG0[j];
	}

//Sコイルの重みを変更

//Gdiff[1][0] *=0.4; 
//Gdiff[3][0] *=0.57; 




/////// 5x1のGdiffを作成
//printf("Gdiff:%lf %lf %lf %lf %lf\n",Gdiff[0][0],Gdiff[1][0],Gdiff[2][0],Gdiff[3][0],Gdiff[4][0]);
////=====行列掛け算  dd=InvCxB * Gdiff 5x5 ==========================
for(i=0;i<5;i++){
dd[i]=0;
	for(j=0;j<1;j++){
//		for(k=0;k<5;k++){
		for (k = 0; k < NUM_SENT_COIL; k++) {
			dd[i]+=InvCxB[i][k]*Gdiff[k][j];
		}
	}
}

 for(i=3;i<5;i++){
 	if(dd[i]> 1000000.0)return(1);
 	if(dd[i]< -1000000.0)return(1);
 }

//printf("dd:%lf %lf %lf %lf %lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
//============ddは1x5横ベクトルにしてる=========================================
//更新ベクトルを算出 delta d　デスよ k:5x5 A':5x6 

//	printf("dd[3]:%lf\n",dd[3]);
angle_wrap1(&dd[3]);       //更新ベクトルのΔθ1、Δθ2の値を丸める
//	printf("dd[4]:%lf\n",dd[4]);
angle_wrap1(&dd[4]);

for(i=0;i<Mlen;i++){ //ここで重みごとに予測位置の更新
	for(j=0;j<3;j++){
    	POSbuff[i][j] =  POS[a][j]+ dd[j] * M[i];     //各重みに対して更新を実行 Dとかddは縦ベクトル　列は更新係数に対応
		}
	}

for(i=0;i<Mlen;i++){ //ここで重みごとに予測位置の更新
	for(j=0;j<2;j++){
    	anglebuff[i][j] =  angles[a][j]+ dd[j+3] * M[i];     //各重みに対して更新を実行 Dとかddは縦ベクトル　列は更新係数に対応
		}
	}
for(i=0;i<Mlen;i++){
	angle_wrap1(&anglebuff[i][0]);
	angle_wrap1(&anglebuff[i][1]);
	}

for(i=0;i<Mlen;i++){
    CC = cos(anglebuff[i][0])*cos(anglebuff[i][1]);
    CS = cos(anglebuff[i][0])*sin(anglebuff[i][1]);
    Si = sin(anglebuff[i][0]);  //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある

	for(j=0;j<NUM_SENT_COIL;j++){
        Dfunc(&POSbuff[i][0],Qm[j],&Scoil_pole[j][0],CC,CS,Si,&estG0[j]);
        //コイルの位置情報を受け取って、その場所での受信信号のゲインの理論値(横ベクトル)を返す関数

		}

	S[i]=0;

//S[i]=pow(Zy[RLPi][0]-estG0[0],2)+0.4 *pow(Zy[RLPi][1]-estG0[1],2)+pow(Zy[RLPi][2]-estG0[2],2)+pow(Zy[RLPi][3]-estG0[3],2)+pow(Zy[RLPi][4]-estG0[4],2)+pow(Zy[RLPi][5]-estG0[5],2);

	


 for(j=0;j<NUM_SENT_COIL;j++){
 		S[i]+=pow(Zy[RLPi][j]-estG0[j],2);
//		S[i]+=fabs(Zy[RLPi][j]-estG0[j]);
		}//各重みで更新したDnewからそれぞれ予測誤差規範を計算


//printf("S[%d]:%lf\n",i,S[i]);
	}


I=min_index(S,0,Mlen);    //予測誤差規範Sが最小になるインデックス番号Iを求める

for(j=0;j<3;j++)POS[a+1][j]=POSbuff[I][j];
for(j=0;j<2;j++)angles[a+1][j]=anglebuff[I][j];
*S_return=S[I];
//予測誤差規範の値


// printf("POS[a+1][j]:%lf %lf %lf \n",POS[a+1][0],POS[a+1][1],POS[a+1][2]);
// printf("angles[a+1][j]:%lf %lf\n",angles[a+1][0],angles[a+1][1]);
// printf("採用するM:%d S[%d]:%lf\n",I,I,S[I]);
// printf("予測誤差規範の値,Calc_ddから返す前:%lf\n",*S_return);
for(i=0;i<NUM_SENT_COIL;i++)
free(A[i]);
free(A);
for(i=0;i<5;i++)
free(B[i]);
free(B);
for(i=0;i<5;i++)
free(C[i]);
free(C);
for(i=0;i<5;i++)
free(InvC[i]);
free(InvC);
for(i=0;i<5;i++)
free(InvCxB[i]);
free(InvCxB);
for(i=0;i<NUM_SENT_COIL;i++)
free(Gdiff[i]);
free(Gdiff);
for(i=0;i<Mlen;i++)
free(POSbuff[i]);
free(POSbuff);
for(i=0;i<Mlen;i++)
free(anglebuff[i]);
free(anglebuff);

return(0);

// (更新された受信コイル推定位置,同角度, 予測誤差規範)
// このプログラム全体を繰り返すことで更新していく
}

int Calc_dd_loop(double **Zy,int RLPi,double **POS,double **angles,int update_max, int *update_num, double **Scoil_pole,double *Qm,double *S_return){
// 更新ベクトル Δdを計算してそれを加え、更新した受信コイル推定位置と角度を返す。
// また、そのときの予測誤差規範Sも返す
// 引数は実測ゲイン:[[0]]xNUM_SENT_COIL  受信コイル位置、角度、送信コイルの位置角度と磁荷NUM_SENT_COILx6
int i,j,k,p,I, i_loop;
double Q[NUM_SENT_COIL];
double estG0[NUM_SENT_COIL],dd[5];
double CC,CS,Si,SC,SS,Ci;

double Rp[3],Rn[3],lenRp,lenRn;
double roRp3x1,roRp3x2,roRp3x3,roRn3x1,roRn3x2,roRn3x3;
double **A,**B,**C,**InvC,**InvCxB;
double **Gdiff;
double M[]={1.0,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005, 0.00001, 5.0e-6, 1.0e-6, 5.0e-7, 1.0e-7, 5.0e-8, 1.0e-8};//更新ベクトルにかける重み
int Mlen =17;
// double M[]={1};
// int Mlen=1;
double **POSbuff,**anglebuff;
double S[17];




A=new_matrix(NUM_SENT_COIL,5);
B=new_matrix(5,NUM_SENT_COIL);
C=new_matrix(5,5);
InvC=new_matrix(SN - 1, SN - 1);
//InvCxB=new_matrix(5,5);
InvCxB = new_matrix(5, NUM_SENT_COIL);
Gdiff=new_matrix(NUM_SENT_COIL,1);
POSbuff=new_matrix(Mlen,3);
anglebuff=new_matrix(Mlen,2);

for(j=0;j<NUM_SENT_COIL;j++){
//printf("%d,Q[j]%lf \n",j,Qm[j]);
	Q[j]=Qm[j]/(4.0*M_PI*(1.26/pow(10.0,4.0)));
//	Q[j]=Q[j]/(16.0*M_PI*M_PI/pow(10,10));//理論どおりではたぶん10
}

for (p = 0; p < update_max; p++)
	S_return[p] = DBL_MAX;//二乗予測誤差規範 初期化

//-----------------------ループ開始 --------------------//

for (i_loop = 1; i_loop < update_max; i_loop++) {
//-----------------------ヤコビ行列算出-------by Okamura--------------------//
////---まず今の受信コイル位置での予測信号ゲインを計算する それを位置情報で偏微分する
CC = cos(angles[i_loop - 1][0])*cos(angles[i_loop - 1][1]);
CS = cos(angles[i_loop - 1][0])*sin(angles[i_loop - 1][1]);
Si = sin(angles[i_loop - 1][0]);  //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある
SC = sin(angles[i_loop - 1][0])*cos(angles[i_loop - 1][1]);
SS = sin(angles[i_loop - 1][0])*sin(angles[i_loop - 1][1]);
Ci = cos(angles[i_loop - 1][0]);


for(i=0;i<NUM_SENT_COIL;i++){ //i番目の送信コイルに対して
//printf("Q[j]%lf \n",Q[j]);	
	for(j=0;j<3;j++){
	Rp[j]=0.0;
	Rn[j]=0.0;
	Rp[j] = POS[i_loop - 1][j] - Scoil_pole[i][j];//i番目送信コイル正極から受信コイル位置までのベクトル
    Rn[j] = POS[i_loop - 1][j] - Scoil_pole[i][j+3];//
	}

lenRp=0.0;
lenRn=0.0;

	for(j=0;j<3;j++){
	lenRp=lenRp+pow(Rp[j],2.0);
	lenRn=lenRn+pow(Rn[j],2.0);   //ベクトルの長さ
	}

	lenRp=pow(lenRp,0.5);
	lenRn=pow(lenRn,0.5);


    roRp3x1 = 3.0 * lenRp * Rp[0];  //偏微分らしいけどなんでこうなるの？やまさんの論文9ページ
    roRp3x2 = 3.0 * lenRp * Rp[1];
    roRp3x3 = 3.0 * lenRp * Rp[2];
    roRn3x1 = 3.0 * lenRn * Rn[0];
    roRn3x2 = 3.0 * lenRn * Rn[1];
    roRn3x3 = 3.0 * lenRn * Rn[2];

  A[i][0] = Q[i] * (  ( (pow(lenRp,3) - Rp[0]*roRp3x1) / (pow(lenRp,6)) - (pow(lenRn,3) - Rn[0]*roRn3x1) / (pow(lenRn,6)) ) * CC           + ( ((-Rp[1])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x1) / (pow(lenRn,6)) ) * CS
                  - ( ((-Rp[2])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x1) / (pow(lenRn,6)) ) * Si  );


  A[i][1] = Q[i] * (  ( ((-Rp[0])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x2) / (pow(lenRn,6)) ) * CC
           +( (pow(lenRp,3) - (Rp[1]*roRp3x2)) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[1])*roRn3x2) / (pow(lenRn,6)) ) * CS
                 - ( ((-Rp[2])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x2) / (pow(lenRn,6)) ) * Si);
             
  A[i][2] = Q[i] * (  ( ((-Rp[0])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x3) / (pow(lenRn,6)) ) * CC
           +( ((-Rp[1])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x3) / (pow(lenRn,6)) ) * CS
                - ( (pow(lenRp,3) - (Rp[2])*roRp3x3) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[2])*roRn3x3) / (pow(lenRn,6)) ) * Si  );

  A[i][3] = Q[i] * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*SC
                 -( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*SS
                  -( (Rp[2])/(pow(lenRp,3)) - (Rn[2])/(pow(lenRn,3)) )*Ci   );
                  
  A[i][4] = Q[i] * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*CS
                 +( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*CC   );  

 //Aはヤコビアン行列（6x5、各行がコイル番号）
//printf("Q[i]%lf \n",Q[i]);

//printf("%lf %lf %lf %lf %lf\n",A[i][0],A[i][1],A[i][2],A[i][3],A[i][4]);
}
//--------------------------------------------------------------------//
//---------------------------受信信号予測値算出-------------------------//



for(i=0;i<NUM_SENT_COIL;i++){
Dfunc(&POS[i_loop - 1][0],Qm[i],&Scoil_pole[i][0],CC,CS,Si,&estG0[i])  ;  //受信信号予測値（NUM_SENT_COILx[[0]]、各行が送信コイル番号に対応）;
}
//--------------------------------------------------------------------//

//printf("B[i][j] \n");
////----転地行列-Bの計算---------------
for(i=0;i<5;i++){
	for(j=0;j<NUM_SENT_COIL;j++){
		B[i][j]=A[j][i];
	}
//printf("%lf %lf %lf %lf %lf %lf\n",B[i][0],B[i][1],B[i][2],B[i][3],B[i][4],B[i][5]);
}
//----------B:5x6-------------------

for (i = 0; i < 5; i++)
	for (j = 0; j < 5; j++)
		C[i][j] = 0.0;

////=====行列掛け算  C 5x5 ========
//printf("C[i][j] \n");
for(i=0;i<5;i++){
	for(j=0;j<5;j++){
		for(k=0;k<NUM_SENT_COIL;k++){
			C[i][j]+=B[i][k]*A[k][j];
		}
	}

//printf("%lf %lf %lf %lf %lf\n",C[i][0],C[i][1],C[i][2],C[i][3],C[i][4]);
}//==============================



//////--逆行列の計算------------------------
Calc_Inv_main(C, InvC);
//InvC=Calc_Inv_main(C); // Calc_Inv_mainは塚本氏によるLU分解を用いた逆行列計算、2は吐き出し法を用いた逆行列計算
//if(InvC==NULL){
//printf("InvC=NULL!\n");
////return((double)DBL_MAX);
//return (0);
//}
//////------------------------------------
//puts("unnko4");

//for (i = 0; i < 5; i++)
//	for (j = 0; j < 5; j++)
//		InvCxB[i][j] = 0.0;
for (i = 0; i < 5; i++)
	for (j = 0; j < NUM_SENT_COIL; j++)
		InvCxB[i][j] = 0.0;

////=====行列掛け算  InvCxB 5x5 ==========================
//for(i=0;i<5;i++){
//	for(j=0;j<5;j++){
//		for(k=0;k<5;k++){
//			InvCxB[i][j]+=InvC[i][k]*B[k][j];
//		}
//	}
//}
for(i=0;i<5;i++){
	for(j=0;j<NUM_SENT_COIL;j++){
		for(k=0;k<5;k++){
			InvCxB[i][j]+=InvC[i][k]*B[k][j];
		}
	}
}
//=====================================================

for(j=0;j<NUM_SENT_COIL;j++){
	Gdiff[j][0]=Zy[RLPi][j]-estG0[j];
	}

//Sコイルの重みを変更

//Gdiff[1][0] *=0.4; 
//Gdiff[3][0] *=0.57; 




/////// 5x1のGdiffを作成
//printf("Gdiff:%lf %lf %lf %lf %lf\n",Gdiff[0][0],Gdiff[1][0],Gdiff[2][0],Gdiff[3][0],Gdiff[4][0]);
////=====行列掛け算  dd=InvCxB * Gdiff 5x5 ==========================
for(i=0;i<5;i++){
dd[i]=0;
	for(j=0;j<1;j++){
//		for(k=0;k<5;k++){
		for (k = 0; k < NUM_SENT_COIL; k++) {
			dd[i]+=InvCxB[i][k]*Gdiff[k][j];
		}
	}
}

 for(i=3;i<5;i++){
// 	if(dd[i]> 1000000.0)return(1);
// if(dd[i]< -1000000.0)return(1);
	 if(dd[i]> 1000000.0)return(0);
	 if(dd[i]< -1000000.0)return(0);
 }

//  for(i=3;i<5;i++){
// 	if(dd[i]> 1000000.0)dd[i] = 0.0;
// 	if(dd[i]< -1000000.0)dd[i] = 0.0;
// }

//printf("dd:%lf %lf %lf %lf %lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
//============ddは1x5横ベクトルにしてる=========================================
//更新ベクトルを算出 delta d　デスよ k:5x5 A':5x6 

//	printf("dd[3]:%lf\n",dd[3]);
angle_wrap1(&dd[3]);       //更新ベクトルのΔθ1、Δθ2の値を丸める
//	printf("dd[4]:%lf\n",dd[4]);
angle_wrap1(&dd[4]);

for(i=0;i<Mlen;i++){ //ここで重みごとに予測位置の更新
	for(j=0;j<3;j++){
    	POSbuff[i][j] =  POS[i_loop - 1][j]+ dd[j] * M[i];     //各重みに対して更新を実行 Dとかddは縦ベクトル　列は更新係数に対応
		}
	}

for(i=0;i<Mlen;i++){ //ここで重みごとに予測位置の更新
	for(j=0;j<2;j++){
    	anglebuff[i][j] =  angles[i_loop - 1][j]+ dd[j+3] * M[i];     //各重みに対して更新を実行 Dとかddは縦ベクトル　列は更新係数に対応
		}
	}
for(i=0;i<Mlen;i++){
	angle_wrap1(&anglebuff[i][0]);
	angle_wrap1(&anglebuff[i][1]);
	}

for(i=0;i<Mlen;i++){
    CC = cos(anglebuff[i][0])*cos(anglebuff[i][1]);
    CS = cos(anglebuff[i][0])*sin(anglebuff[i][1]);
    Si = sin(anglebuff[i][0]);  //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある

	for(j=0;j<NUM_SENT_COIL;j++){
        Dfunc(&POSbuff[i][0],Qm[j],&Scoil_pole[j][0],CC,CS,Si,&estG0[j]);
        //コイルの位置情報を受け取って、その場所での受信信号のゲインの理論値(横ベクトル)を返す関数

		}

	S[i]=0;

//S[i]=pow(Zy[RLPi][0]-estG0[0],2)+0.4 *pow(Zy[RLPi][1]-estG0[1],2)+pow(Zy[RLPi][2]-estG0[2],2)+pow(Zy[RLPi][3]-estG0[3],2)+pow(Zy[RLPi][4]-estG0[4],2)+pow(Zy[RLPi][5]-estG0[5],2);

	


 for(j=0;j<NUM_SENT_COIL;j++){
 		S[i]+=pow(Zy[RLPi][j]-estG0[j],2);
//		S[i]+=fabs(Zy[RLPi][j]-estG0[j]);
		}//各重みで更新したDnewからそれぞれ予測誤差規範を計算


//printf("S[%d]:%lf\n",i,S[i]);
	}


I=min_index(S,0,Mlen);    //予測誤差規範Sが最小になるインデックス番号Iを求める

for(j=0;j<3;j++)POS[i_loop][j]=POSbuff[I][j];
for(j=0;j<2;j++)angles[i_loop][j]=anglebuff[I][j];
S_return[i_loop]=S[I];//予測誤差規範の値

if(S_return[i_loop] >= S_return[i_loop - 1])
	break;

}

*update_num = i_loop;


// printf("POS[a+1][j]:%lf %lf %lf \n",POS[a+1][0],POS[a+1][1],POS[a+1][2]);
// printf("angles[a+1][j]:%lf %lf\n",angles[a+1][0],angles[a+1][1]);
// printf("採用するM:%d S[%d]:%lf\n",I,I,S[I]);
// printf("予測誤差規範の値,Calc_ddから返す前:%lf\n",*S_return);
for(i=0;i<NUM_SENT_COIL;i++)
free(A[i]);
free(A);
for(i=0;i<5;i++)
free(B[i]);
free(B);
for(i=0;i<5;i++)
free(C[i]);
free(C);
for(i=0;i<5;i++)
free(InvC[i]);
free(InvC);
for(i=0;i<5;i++)
free(InvCxB[i]);
free(InvCxB);
for(i=0;i<NUM_SENT_COIL;i++)
free(Gdiff[i]);
free(Gdiff);
for(i=0;i<Mlen;i++)
free(POSbuff[i]);
free(POSbuff);
for(i=0;i<Mlen;i++)
free(anglebuff[i]);
free(anglebuff);

return(0);

// (更新された受信コイル推定位置,同角度, 予測誤差規範)
// このプログラム全体を繰り返すことで更新していく
}
void ten0(int P,double **tens){

 int PP,PP2,PP3;
 int i,j,k,a;
 double LL;

	if(P==5)LL = 100.0;  //測定領域の一辺の長さ[m]
	if(P==3)LL = 100.0;  //測定領域の一辺の長さ[m]
	if(P==7 || P==9)LL = 100.0;//測定領域の一辺の長さ[m]
	//if(P==75 || P==45){L = 80.0;L2 = 100.0;L3 = 100.0;}

	if(P==5||P==3){
		PP = (P-1)/2;
		a = 0;
		for(i=-PP;i<=PP;i++){
			for(j=-PP;j<=PP;j++){
				for(k=-PP;k<=PP;k++){
					tens[a][0] = (double)k / PP;
					tens[a][1] = (double)j / PP;
					tens[a][2] = (double)i / PP;
					a +=1;
				}
			}
		}
	
		for(i=0;i<(int)pow((double)P,3.0);i++){
			for(j=0;j<3;j++){
				tens[i][j] = tens[i][j] * LL / 2.0;	
			}
		}
	}else if(P==7){
	
		tens[0][0]=0;tens[0][1]=0;tens[0][2]=-1;
		tens[1][0]=0;tens[1][1]=-1;tens[1][2]=0;
		tens[2][0]=-1;tens[2][1]=0;tens[2][2]=0;
		tens[3][0]=0;tens[3][1]=0;tens[3][2]=0;
		tens[4][0]=1;tens[4][1]=0;tens[4][2]=0;
		tens[5][0]=0;tens[5][1]=1;tens[5][2]=0;
		tens[6][0]=0;tens[6][1]=0;tens[6][2]=1;
	
		for(i=0;i<P;i++){
			for(j=0;j<3;j++){
				tens[i][j] = tens[i][j] * LL / 2.0;
			}
		}
	}else if(P==9){
	
		tens[0][0]=-1;tens[0][1]=-1;tens[0][2]=-1;
		tens[1][0]=1;tens[1][1]=-1;tens[1][2]=-1;
		tens[2][0]=-1;tens[2][1]=1;tens[2][2]=-1;
		tens[3][0]=1;tens[3][1]=1;tens[3][2]=-1;
		tens[4][0]=0;tens[4][1]=0;tens[4][2]=0;
		tens[5][0]=-1;tens[5][1]=-1;tens[5][2]=1;
		tens[6][0]=1;tens[6][1]=-1;tens[6][2]=1;
		tens[7][0]=-1;tens[7][1]=1;tens[7][2]=1;
		tens[8][0]=1;tens[8][1]=1;tens[8][2]=1;
	
		for(i=0;i<P;i++){
			for(j=0;j<3;j++){
				tens[i][j] = tens[i][j] * LL / 2.0;
			}
		}
 
		//for(a = 0; a < P; a ++){
		//	rotation(tens[a][0],tens[a][1],tens[a][2],&mx,&my,&mz);
		//	hlmmmtocount(tens[a][0],tens[a][1],tens[a][2],&mx,&my,&mz);
		//	tens[a][0] = mx;
		//	tens[a][1] = my;
		//	tens[a][2] = mz;
		//	}

 
	}else{
//		printf("渡されてるPの値がおかしい @ten0.c P:%d\n",P);
//		exit(1);
		for (i = 0; i < P; i++)
			for (j = 0; j < 3; j++)
				tens[i][j] = 0.0;
	}
/*
//-----確認用----
//printf("CP:%d\n",P);
// printf("L:%lf a:%d\n",L,a);
for(i=0;i<a;i++){
	//printf("%lf %lf %lf\n", tens[i][0],tens[i][1],tens[i][2]);
 }
// //[m]で返す
// printf("unnko\n");
////------------
*/
}

void ten189(int P, int i_ch, double **tens){
	int a, i, i_x, i_y, i_z, j;

	a = 0;

	if (P == 189)
		for (i_z = -1; i_z <= 1; i_z++)
			for (i_y = -4; i_y <= 4; i_y++)
				for (i_x = -3; i_x <= 3; i_x++) {
					tens[a][0] = (i_x + (i_ch - 1) / 4 - 1) * 25.0;
					tens[a][1] = (i_y + (i_ch - 1) % 4 - 1) * 25.0;
					tens[a][2] = i_z * 50.0;
					a++;
				}
	else {
		//		printf("渡されてるPの値がおかしい @ten0.c P:%d\n",P);
		//		exit(1);
		for (i = 0; i < P; i++)
			for (j = 0; j < 3; j++)
				tens[i][j] = 0.0;
	}
	/*
	//-----確認用----
	//printf("CP:%d\n",P);
	// printf("L:%lf a:%d\n",L,a);
	for(i=0;i<a;i++){
	//printf("%lf %lf %lf\n", tens[i][0],tens[i][1],tens[i][2]);
	}
	// //[m]で返す
	// printf("unnko\n");
	////------------
	*/
}

void SDM_2(double **POS,double *dd,int a,double *Scoil_pole,double *A){

int i,j,k,I;
double Q;
double estG0[NUM_SENT_COIL];
double CC,CS,Si,SC,SS,Ci;

double Rp[3],Rn[3],lenRp,lenRn;
double roRp3x1,roRp3x2,roRp3x3,roRn3x1,roRn3x2,roRn3x3;
double **POSbuff,**anglebuff;




    Q = 1.0/(4.0*M_PI*(1.26/pow(10.0,4.0)));





CC = cos(dd[3]/180.0*M_PI)*cos(dd[4]/180.0*M_PI);
CS = cos(dd[3]/180.0*M_PI)*sin(dd[4]/180.0*M_PI);
Si = sin(dd[3]/180.0*M_PI);  //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある
SC = sin(dd[3]/180.0*M_PI)*cos(dd[4]/180.0*M_PI);
SS = sin(dd[3]/180.0*M_PI)*sin(dd[4]/180.0*M_PI);
Ci = cos(dd[3]/180.0*M_PI);


    for(j=0;j<3;j++){
    Rp[j]=0.0;
    Rn[j]=0.0;
    Rp[j] = POS[a][j]+dd[j] - Scoil_pole[j];//i番目送信コイル正極から受信コイル位置までのベクトル
    Rn[j] = POS[a][j]+dd[j] - Scoil_pole[j+3];//
    }

lenRp=0.0;
lenRn=0.0;

    for(j=0;j<3;j++){
    lenRp=lenRp+pow(Rp[j],2.0);
    lenRn=lenRn+pow(Rn[j],2.0);   //ベクトルの長さ
    }

    lenRp=pow(lenRp,0.5);
    lenRn=pow(lenRn,0.5);


    roRp3x1 = 3.0 * lenRp * Rp[0];  //偏微分らしいけどなんでこうなるの？やまさんの論文9ページ
    roRp3x2 = 3.0 * lenRp * Rp[1];
    roRp3x3 = 3.0 * lenRp * Rp[2];
    roRn3x1 = 3.0 * lenRn * Rn[0];
    roRn3x2 = 3.0 * lenRn * Rn[1];
    roRn3x3 = 3.0 * lenRn * Rn[2];

  A[0] = Q * (  ( (pow(lenRp,3) - Rp[0]*roRp3x1) / (pow(lenRp,6)) - (pow(lenRn,3) - Rn[0]*roRn3x1) / (pow(lenRn,6)) ) * CC
			+ ( ((-Rp[1])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x1) / (pow(lenRn,6)) ) * CS
                  - ( ((-Rp[2])*roRp3x1) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x1) / (pow(lenRn,6)) ) * Si  );


  A[1] = Q * (  ( ((-Rp[0])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x2) / (pow(lenRn,6)) ) * CC
           +( (pow(lenRp,3) - (Rp[1]*roRp3x2)) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[1])*roRn3x2) / (pow(lenRn,6)) ) * CS
                 - ( ((-Rp[2])*roRp3x2) / (pow(lenRp,6)) - ((-Rn[2])*roRn3x2) / (pow(lenRn,6)) ) * Si);
             
  A[2] = Q * (  ( ((-Rp[0])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[0])*roRn3x3) / (pow(lenRn,6)) ) * CC
           +( ((-Rp[1])*roRp3x3) / (pow(lenRp,6)) - ((-Rn[1])*roRn3x3) / (pow(lenRn,6)) ) * CS
                - ( (pow(lenRp,3) - (Rp[2])*roRp3x3) / (pow(lenRp,6)) - (pow(lenRn,3) - (Rn[2])*roRn3x3) / (pow(lenRn,6)) ) * Si  );

  A[3] = Q * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*SC
                 -( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*SS
                  -( (Rp[2])/(pow(lenRp,3)) - (Rn[2])/(pow(lenRn,3)) )*Ci   );
                  
  A[4] = Q * ( -( (Rp[0])/(pow(lenRp,3)) - (Rn[0])/(pow(lenRn,3)) )*CS
                 +( (Rp[1])/(pow(lenRp,3)) - (Rn[1])/(pow(lenRn,3)) )*CC   ); 

  //printf("A %lf %lf %lf %lf %lf\n",A[0],A[1],A[2],A[3],A[4]);

}



int main(void)
{
	int i, i_ch, i_a, i_p, j, k, a, p, b, old_b;
	int zure_suitei, zerocorrection;
	int file_no, i_data_no;
	int evaluation = 1; // 0:ただの位置推定 1:推定精度評価 Scoilfilename GPfilenameを使ってDatafilenameのを評価するだからDatafilenameは125点のデータね
	int Lines, DP;
	int I2;
	int update_max1 = 50, update_max = 10000;
	int IAlen = (int)pow(17.0,2.0);
	int IA;
	int PP, PP2, l, U_N;
	int ss;
	int LOOP, MIN = 0;
	int MAX;
	int m;
	int u;
	int aiueo;
	int *L_error_index;
	int update_num[PATERN];
	char Sc_dir[FILENAME_MAX];
	char Data_dir[FILENAME_MAX];
	char Datafilename[FILENAME_MAX];
	char ScoilDirname[FILENAME_MAX] = ""; //"ga9kai_0330_new_calib";
	char Scoilfilename[FILENAME_MAX]; //"gakkai_0330_new_calib";
	char GPfilename[FILENAME_MAX]; //"gakkai_0330_new_calib";
	char Spfilename[FILENAME_MAX];
	char resultflname[FILENAME_MAX];
	char out_flname[FILENAME_MAX];
	char out_fldir[FILENAME_MAX];
	double check_state;
	double mean_signal_error;
	double POSerrormean;
	double s0, sum;
	double res_error_mean;
	double sum_buf_10000;
	double r[NUM_SENT_COIL];
	double res[3];
	double alpha[PATERN];
	double sum2[PATERN], sum3[PATERN];
	double d_f[5];
	double angle_av[2] = {0.0,0.0};
	double dd_buf[PATERN][5];
	double test_box[TEST][5];
	double A2[NUM_SENT_COIL][189][5];
	double *Qm;
	double *signal_error, *signal_error2, *mean_signal_error2, *S, *S1;
	double *POSerror;
	double *res_error;
	double *dd, *dd2;
	double *max_buf;
	double **dd_stock;
	double **Scoil;
	double **dd_read;
	double **data_buf;
	double **RLP;
	double **pos, **angles, **POSAmS, **Scoil_pole;
	double **Initial_Angle;
	double **U_B;
	double **Z;
	double **dummy;
	FILE *fp;
	FILE *rslt_fp;

	sprintf(Sc_dir, "%s", SCOIL_DIR);

	for (i_data_no = START_DATA_NO; i_data_no <= END_DATA_NO; i_data_no++) {
		sprintf(Data_dir, "c-mean");

		//	printf("put Scoil file dir\n");
		//	scanf("%s", Sc_dir);
		//	printf("put Data file dir\n");
		//	scanf("%s", Data_dir);
		//	printf("1.Set_starting_point?      ->enter:1\n2.Read_location_file       ->enter:2\n3.Normal_estimation        ->enter:3\n");
		//	scanf("%d",&aiueo);
		aiueo = 3;
		//	printf("zure_suitei?(O.K. = 1, No = other)\n");
		//	scanf("%d", &zure_suitei);
		zure_suitei = 0;
		//	printf("zero_correction?(O.k. = 1, No = other)\n");
		//	scanf("%d", &zerocorrection);
		zerocorrection = 0;

		for (i_ch = START_CH_NO; i_ch <= END_CH_NO; i_ch++) {
			if (i_ch == REF_CH_NO)
				continue;

			U_B = new_matrix(9 * 25, 5);
			PP2 = 2;
			PP = 2;

			sprintf(resultflname, "..\\..\\data\\%s\\%s\\result_file_ch%d.data", DATA_DATE, Data_dir, i_ch);
			if ((rslt_fp = fopen(resultflname, "a")) == NULL) {
				printf("resultflname = %s\n", resultflname);
				printf("file open error @resultflname\n");
				return EXIT_FAILURE;
			}
			//puts("ok");//////////////////////////////////////////////////////////
			////U_Bは位置推定の初期位置を決める行列＊今は使ってない////////
			//////////////////////////////////////////////////////////// /
			for (k = 0; k < 9; k++) {
				a = 0;
				for (i = -PP2; i <= PP2; i++)
					for (j = -PP; j <= PP; j++) {
						U_B[a + (k * (PP * 2 + 1) * (PP2 * 2 + 1))][3] = ((double)i / PP2) * 20.0;
						U_B[a + (k * (PP * 2 + 1) * (PP2 * 2 + 1))][4] = ((double)j / PP) *20.0;
						a += 1;
					}
			}
			a = 0;
			for (m = -4; m <= 4; m++)
				for (l = -4; l <= 4; l++)
					for (k = -1; k <= 1; k++)
						for (j = -1; j <= 1; j++)
							for (i = -1; i <= 1; i++) {
								test_box[a][0] = (double)i * 40.0;
								test_box[a][1] = (double)j * 40.0;
								test_box[a][2] = (double)k * 40.0;
								test_box[a][3] = ((double)l / 4.0) * 40.0;
								test_box[a][4] = ((double)m / 4.0) * 40.0;
								a++;
								//printf("%lf\t%lf\n",((double)l/4.0),((double)m/4.0));
							}
			l = 0;
			for (i = -1; i <= 1; i += 2)
				for (j = -1; j <= 1; j += 2)
					for (k = -1; k <= 1; k += 2) {
						for (a = 0; a < (PP * 2 + 1) * (PP2 * 2 + 1); a++) {
							if (l == 0) {
								U_B[a][0] = 0.0;
								U_B[a][1] = 0.0;
								U_B[a][2] = 0.0;
							}
							//	printf("%d\n",a+l*(PP*2+1)*(PP2*2+1));
							U_B[a + (l + 1) * (PP * 2 + 1) * (PP2 * 2 + 1)][0] = (double)i * 20.0;
							U_B[a + (l + 1) * (PP * 2 + 1) * (PP2 * 2 + 1)][1] = (double)j * 20.0;
							U_B[a + (l + 1) * (PP * 2 + 1) * (PP2 * 2 + 1)][2] = (double)k * 20.0;
						}
						l++;
					}
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//////////////////////////////////////////////////////////////////////////////////////////////////////////	
			puts("ok");
			//		scanf("%d", &check_state);
					//// memory allication
			Scoil = new_matrix(NUM_SENT_COIL, 6);
			dd_read = new_matrix(FILENUM, AXIS_NUM);
			Qm = (double *)malloc(sizeof(double) * NUM_SENT_COIL);
			for (file_no = 0; file_no < DATA_FILENUM; file_no++) {
				sprintf(Datafilename, "..\\..\\data\\%s\\%s\\out_ch%d", DATA_DATE, Data_dir, i_ch);
				printf("%s\n", Datafilename);
				sprintf(Scoilfilename, "..\\..\\data\\%s\\c-mean\\out_ch%d", Sc_dir, i_ch);
				sprintf(GPfilename, "..\\..\\data\\%s\\c-mean\\out_ch%d", Sc_dir, i_ch);
				sprintf(Spfilename, "..\\..\\data\\%s\\c-mean\\out_ch%d", Sc_dir, i_ch);
				printf("test6\n");
				data_buf = Data_import(Datafilename, ScoilDirname, Scoilfilename, GPfilename, Spfilename, &Lines, Scoil, Qm, r, dd_read);//データセット・送信コイル座標・ゲインファイルのオープン
				printf("test5\n");
				printf("%d:Lines\n", Lines);
				////// 確認用--------------------	
				for (i = 0; i < Lines; i++) {
					printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", data_buf[i][0], data_buf[i][1], data_buf[i][2], data_buf[i][3], data_buf[i][4], data_buf[i][5], data_buf[i][6], data_buf[i][7]);
				}
				printf("%d:Lines\n", Lines);
				puts("Scoil_data");
				for (i = 0; i < NUM_SENT_COIL; i++) {
					printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Scoil[i][0], Scoil[i][1], Scoil[i][2], Scoil[i][3], Scoil[i][4], Scoil[i][5], Qm[i]);
				}
				puts("\nScoil_data");
				puts("");
				/////----------------------------
				for (i = 0; i < Lines; i++)
					for (j = 0; j < NUM_SENT_COIL; j++)
						data_buf[i][j] -= r[j];
				if (evaluation == 1)
					Lines = Lines - 1;
				//////---memory allocation------------------
				pos = new_matrix(update_max, 3);
				angles = new_matrix(update_max, 2);
				signal_error = (double *)malloc(sizeof(double) * Lines);
				signal_error2 = (double *)malloc(sizeof(double) * Lines);
				mean_signal_error2 = (double *)malloc(sizeof(double) * IAlen);
				POSAmS = new_matrix(Lines, 5);
				S = (double *)malloc(sizeof(double) * update_max);
				S1 = (double *)malloc(sizeof(double) * IAlen);
				L_error_index = (int *)malloc(sizeof(int) * Lines);
				//////----------------------------------------
				/////--送信コイルの情報----------------------
				//送信コイルの極の座標を返す関数 6x6
				Scoil_pole = new_matrix(NUM_SENT_COIL, 6);
				for (i = 0; i < NUM_SENT_COIL; i++)
					for (j = 0; j < 6; j++)
						Scoil_pole[i][j] = Scoil[i][j];
				for (i = 0; i < NUM_SENT_COIL; i++)
					printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Scoil_pole[i][0], Scoil_pole[i][1], Scoil_pole[i][2], Scoil_pole[i][3], Scoil_pole[i][4], Scoil_pole[i][5]);
				////////////////////////-----位置推定の前準備-----受信コイル初期角度の設定--------------------------------
				////初期角度は使ってない＊内田/////////////////////////////
				////初期角度の候補を20degづつ刻む
				Initial_Angle = new_matrix(IAlen, 2);
				for (i = 0; i < pow(IAlen, 0.5); i++)
					for (j = 0; j < pow(IAlen, 0.5); j++) {
						Initial_Angle[i * (int)pow(IAlen, 0.5) + j][0] = (11.0 * j) * M_PI / 180.0;
						Initial_Angle[i * (int)pow(IAlen, 0.5) + j][1] = (11.0 * i) * M_PI / 180.0;
					}//仰角が90°のところは失敗するよ。だから90°にならんように設定せないかん
				////-確認用----
				// 	for(i=0;i<IAlen;i++){
				// 		printf("%dInitial_Angle:%lf %lf\n",i,Initial_Angle[i][0]*180/M_PI,Initial_Angle[i][1]*180/M_PI);
				// 	}
				//////////----------------------------
				////////////////////////-----初期角度の選択と決定---------------------------------
				//初期化--------
				for (i = 0; i < Lines; i++)
					L_error_index[i] = -1;
				sum = 0.0;
				a = 0;
				b = old_b = 0;
				res_error = (double *)malloc(sizeof(double) * Lines);
				max_buf = (double *)malloc(sizeof(double) * Lines);
				MAX = max_index(max_buf, 0, Lines);
				//------------	
				if (evaluation == 1) {
					switch (Lines) {
					case 7: DP = 7;
						break;
					case 9: DP = 9;
						break;
					case 27: DP = 3;
						break;
					case 75: DP = 75;
						break;
					case 125: DP = 5;
						break;
					default: DP = Lines;
						break;
					}
					printf("lins%d\n", Lines);
					POSerror = (double *)malloc(sizeof(double) * Lines);
					RLP = new_matrix(Lines, 3);
					dummy = new_matrix(Lines, 3);
					if (DP == 189)
						ten189(DP, i_ch, RLP);
					else
						ten0(DP, RLP);
					dd = (double *)malloc(sizeof(double) * 5);
					dd2 = (double *)malloc(sizeof(double) * 5);
					dd_stock = new_matrix(PATERN, 5);
					j = 1;
					for (i = 0; i < PATERN; i++) {
						if (i % 4 == 0)
							j++;
						alpha[i] = 1.0 / pow(10.0, j + 3.0) * (double)(4.0 - (i % 4)) * 2.5;
						//printf("%f\n",alpha[i]);
					}
					Z = new_matrix(NUM_SENT_COIL, Lines);
					for (i = 0; i < 5; i++)
						dd[i] = 0.0;
					///ここから最急降下法で受信コイルの原点からのズレを推定/////////////
					if (zure_suitei == 1) {
						for (LOOP = 0; LOOP < 200000; LOOP++) {
							for (i = 0; i < NUM_SENT_COIL; i++) {
								for (k = 0; k < Lines; k++) {
									//dd[0]=0;dd[1]=0;dd[2]=0;
									//printf("DD:%lf\t%lf\t%lf\t%lf\t%lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
									Dfuncg_dd(RLP[k], 1.0, Scoil[i], Z[i], k, i, dd);
									//printf("Z[%d]:%lf\n",k,Z[i][k]);
									SDM_2(RLP, dd, k, Scoil[i], A2[i][k]);
									//printf("Z[%d]:%lf\n",k,Z[i][k]);
									//printf("A %lf %lf %lf %lf %lf\n",A2[i][k][0],A2[i][k][1],A2[i][k][2],A2[i][k][3],A2[i][k][4]);
								}
							}
							//exit(0);
							for (j = 0; j < 5; j++)
								d_f[j] = 0.0;
							for (i = 0; i < NUM_SENT_COIL; i++)
								for (k = 0; k < Lines; k++)
									for (j = 0; j < 5; j++)
										d_f[j] += -2.0 * (data_buf[k][i] - Qm[i] * Z[i][k]) * (Qm[i] * A2[i][k][j]);
							//printf("d_f[%f]\n",d_f[j]);
							//printf("d_buf[%f]",data_buf[k][i]);
							for (k = 0; k < PATERN; k++) {
								for (j = 0; j < 5; j++)
									dd_buf[k][j] = dd[j] - alpha[k] * d_f[j];
								//printf("dd_buf:%lf:%lf:%lf:%lf:%lf\n",dd_buf[k][0],dd_buf[k][1],dd_buf[k][2],dd_buf[k][3],dd_buf[k][4]);
								sum2[k] = 0.0;
								for (i = 0; i < NUM_SENT_COIL; i++)
									for (j = 0; j < Lines; j++)
										Dfuncg_dd(RLP[j], 1.0, Scoil[i], Z[i], j, i, dd_buf[k]);
								//printf("Z[%d]:%lf\n",j,Z[i][j]);
								//printf("dd_buf:%lf:%lf\n",dd_buf[k][3],dd_buf[k][4]);
								//printf("[%d][%d]\t%f\n",k,j,Q[i]*Z[i][j]);
								for (i = 0; i < NUM_SENT_COIL; i++)
									for (j = 0; j < Lines; j++)
										//printf("[%d][%d]\t%f\t%f\n",k,j,calib_buf[i][j],Q[i]*Z[i][j]);
										sum2[k] += pow((data_buf[j][i] - Qm[i] * Z[i][j]), 2.0);
								//sum3[k] += pow((data_buf[j][i]- Qm[i]*Z[i][j])/Qm[i],2.0);
								//sum2[k] += pow((calib_buf[i][j+(N*calib_point)]- Q[i]*Z[i][j+(N*calib_point)]),2.0);
								//sum2[k]=sum2[k]/2.0;
								//printf("sum2[%d]\t%f\n",k,sum2[k]);
							}
							MIN = min_index(sum2, 0, PATERN);
							if (LOOP == 0)
								printf("sum2:%lf\n", sum2[MIN]);
							for (j = 0; j < 5; j++) {
								dd[j] = dd[j] - alpha[MIN] * d_f[j];
							}
							//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
							//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
							if (LOOP == 0)
								sum_buf_10000 = sum2[MIN];
							//printf("%lf\n",sum2[MIN]);
							//printf("%lf\n",sum_buf_10000);
							if (LOOP > 999)
								if (LOOP % 1000 == 0) {
									if ((sum2[MIN] / sum_buf_10000) > 0.9999999999) {
										puts("LOOOPEND\n");
										break;
									}
									//printf("%lf\n",(sum2[MIN]/sum_buf_10000));
									sum_buf_10000 = sum2[MIN];
								}
						}
						printf("d_f:%lf:%lf:%lf:%lf:%lf\n", dd[0], dd[1], dd[2], dd[3], dd[4]);
						printf("sum2:%lf\n", sum2[MIN]);
						//printf("sum2:%lf\n",sum2[MIN]);
						//scanf("%d",&j);
					}
					else {
						for (i = 0; i < 5; i++)
							dd[i] = dd_read[0][i];
					}
					for (i = 0; i < 5; i++)
						printf("dd[%d] = %f\n", i, dd[i]);
					//			scanf("%f", &check_state);
								/////////////////////////////////////////////////////////////////////
								//////////////////////////////////////////////////////////////
					for (i = 0; i < Lines; i++) {//位置推定スタート
						/////////ここは岡村さんの名残・つかってない///////////////////////////////////////////////////////////////////
						//printf("%d:b\n",b);
						if (b != 1 || i == 0) {//b==1は前回が成功のとき だから
							for (j = 0; j < 3; j++) {
								pos[0][j] = 0.0;//一番最初と,前の点で位置推定失敗したとき
							}
							////----初期角度の決定--------------
							for (k = 0; k < IAlen; k++) {
								for (j = 0; j < 2; j++) {
									angles[0][j] = Initial_Angle[k][j]; //初期位置
								}
								//  	angles[0][0]=0;
								//	angles[0][1]=-M_PI/4;
								//		printf("Point:%d Initial_Angle[%d][j]:%lf %lf \n",i,k,angles[0][0]*180/M_PI,angles[0][1]*180/M_PI);
								for (p = 0; p < update_max; p++)
									S[p] = DBL_MAX;//二乗予測誤差規範 初期化
								for (j = 1; j < update_max1; j++) {//j:更新の繰り返し回数
									//				Calc_dd(Zy,i,POS,angles,j-1,Scoil_pole,Qm,&S[j]);
									if (Calc_dd(data_buf, i, pos, angles, j - 1, Scoil_pole, Qm, &S[j]) != 0) {
										//					printf("ddでか杉でbreak!更新%d回目\n",j);
										//					S[j]=DBL_MAX;
										break;
									}
									//				printf("予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[j]);
									if ((j == 2 && S[j] > 1000000000) || (j == 3 && S[j] > 800000000) || (j == 5 && S[j] > 500000000) || (j == 40 && S[j] > 50000000)) {
										printf("break@j=%d point%d\n", j, i);
										break;      //　更新が50回いってもまだでかいようなら次の初期角度での計算に移る
									}
									if (S[j] >= S[j - 1])
										break;
									I2 = min_index(S, 0, update_max1);
									S1[k] = S[I2];//S1はIAlenだけ確保してある
								}
								//printf("S1:%lf %lf %lf %d IA,k\n",S1[k],Initial_Angle[k][0],Initial_Angle[k][1],k);
							}
							/////-----初期角度番号 IA 決定 代入--------
							IA = min_index(S1, 0, IAlen);
							//			for(p=0;p<3;p++)pos[0][p]=0.0; //初期位置
							for (p = 0; p < 2; p++)
								angles[0][p] = Initial_Angle[IA][p]; //初期位置
						}
						else {
							////////////////// b==1 まえの位置推定が成功してるとき
							for (p = 0; p < 3; p++)
								pos[0][p] = POSAmS[i - 1][p]; //初期位置
							for (p = 0; p < 2; p++)
								angles[0][p] = POSAmS[i - 1][p + 3]; //初期位置
						}
						//		printf("Point%d	Initial Angle:%lf[deg] %lf[deg]\n",i,angles[0][0]*180.0/M_PI,angles[0][1]*180.0/M_PI);
						// 		angles[0][0]=0.0*M_PI/180.0;
						//		angles[0][1]=0.0*M_PI/180.0;
						//////----ここまで初期角度の決定--------------
						//SDM_2(tens,dd[N],k,S_p[i],A2[N][i][k+N*calib_point]);
						/*
						for(k=0;k<TEST;k++)
						{
						test_res[k]=0.0;
						test_dd[0]=0.0;
						test_dd[1]=0.0;
						test_dd[2]=0.0;
						test_dd[3]=test_box[k][3];
						test_dd[4]=test_box[k][4];
						for(j=0;j<6;j++){
						Dfuncg_dd(test_box[k],Qm[j],Scoil[j],test_z[j],k,j,test_dd);
						test_res[k] += pow((data_buf[i][j]- test_z[j][k]),2.0);
						}
						}
						for(j=0;j<10;j++){
						m=min_index(test_res,0,TEST);
						for(l=0;l<5;l++){
						//U_B[j][l]=test_box[m][l];
						}
						//printf("test:%lf\t%d\n",test_res[m],m);
						//printf("%lf\t%lf\t%lf\t%lf\t%lf\n",test_box[m][0],test_box[m][1],test_box[m][2],test_box[m][3],test_box[m][4]);
						test_res[m]=DBL_MAX;
						}
						//test_res[k] += pow((data_buf[i][j]- Z[k][j]),2.0);
						*/
						//////////////////////////////////////////////////////////////////////////////////////////////////////
						////////---------位置推定のメイン-----------------------------------
						for (U_N = 0; U_N < 1/*225*/; U_N++) {
							if (U_N != 0) {
								pos[0][0] = U_B[U_N - 1][0];
								pos[0][1] = U_B[U_N - 1][1];//初期位置
								pos[0][2] = U_B[U_N - 1][2];
								for (p = 0; p < 2; p++)
									angles[0][p] = U_B[U_N - 1][p + 3] * M_PI / 180.0;
							}
							//for(p=0;p<3;p++)pos[0][p]=U_B[i][p];
							//angles[0][0]=(-37.971880+180.0)	*M_PI/180.0;
							//angles[0][1]=-41.614744*M_PI/180.0; 
							//updata_max=1;
							if (b != 1)
								for (u = 0; u < 5; u++)
									dd2[u] = 0.0;
							if (b == 3) {
								dd2[0] = 50.0;
								dd2[1] = -50.0;
								dd2[2] = 50.0;
								dd2[3] = -30.0 * M_PI / 180.0;
								dd2[4] = -30.0 * M_PI / 180.0;
							}
							if (b == 4) {
								dd2[0] = 50.0;
								dd2[1] = 50.0;
								dd2[2] = 30.0;
								dd2[3] = 30.0 * M_PI / 180.0;
								dd2[4] = -30.0 * M_PI / 180.0;
							}
							if (b == 5) {
								dd2[0] = -30.0;
								dd2[1] = 30.0;
								dd2[2] = 50.0;
								dd2[3] = 30.0 * M_PI / 180.0;
								dd2[4] = -30.0 * M_PI / 180.0;
							}
							if (b == 6) {
								dd2[0] = -50.0;
								dd2[1] = -50.0;
								dd2[2] = 40.0;
								dd2[3] = -30.0 * M_PI / 180.0;
								dd2[4] = -30.0 * M_PI / 180.0;
							}
							if (b == 7) {
								dd2[0] = 30.0;
								dd2[1] = 50.0;
								dd2[2] = -50.0;
								dd2[4] = -30.0 * M_PI / 180.0;
							}
							if (b == 8) {
								dd2[0] = -10.0;
								dd2[1] = -30.0;
								dd2[2] = -50.0;
								dd2[4] = 30.0 * M_PI / 180.0;
							}
							if (b == 9) {
								dd2[0] = 50.0;
								dd2[1] = -50.0;
								dd2[2] = -50.0;
								dd2[4] = 30.0 * M_PI / 180.0;
							}
							if (b == 10) {
								dd2[0] = -50.0;
								dd2[1] = -10.0;
								dd2[2] = -20.0;
								dd2[3] = 30.0 * M_PI / 180.0;
								dd2[4] = 30.0 * M_PI / 180.0;
							}
							if (b == 11) {
								dd2[0] = -100.0;
								dd2[1] = 125.0;
								dd2[2] = -50.0;
								dd2[3] = 0.0;
								dd2[4] = 0.0;
							}
							if (b == 12) {
								dd2[0] = 0.0;
								dd2[1] = -100.0;
								dd2[2] = 50.0;
								dd2[3] = 0.0;
								dd2[4] = 0.0;
							}
							if (b == 13) {
								dd2[0] = 20.0;
								dd2[1] = 140.0;
								dd2[2] = -160.0;
								dd2[3] = 0.0;
								dd2[4] = 0.0;
							}
							if (b == 14) {
								dd2[0] = 5.0;
								dd2[1] = 160.0;
								dd2[2] = -160.0;
								dd2[3] = -5.0 * M_PI / 180.0;
								dd2[4] = 0.0;
							}
							if (b == 15) {
								dd2[0] = 30.0;
								dd2[1] = 215.0;
								dd2[2] = -260.0;
								dd2[3] = 25.0 * M_PI / 180.0;
								dd2[4] = -15.0 * M_PI / 180.0;
							}
							if (b == 16) {
								dd2[0] = 6.870704;
								dd2[1] = -55.214141;
								dd2[2] = -15.524439;
								dd2[3] = 0.294515;
								dd2[4] = 2.985793;
							}							
							if (aiueo == 1) {
								printf("enter_starting_point:%d\n", i + 1);
								scanf("%lf", &dd2[0]);
								scanf("%lf", &dd2[1]);
								scanf("%lf", &dd2[2]);
							}
							else if (aiueo == 2) {
								double location[27][3] = { -50, -50, -50,
									0, -50, -50,
									50, -50, -50,
									-50, 0, -50,
									0, 0, -50,
									50, 0, -50,
									-50, 50, -50,
									0, 50, -50,
									50, 50, -50,
									-50, -50, 0,
									0, -50, 0,
									50, -50, 0,
									-50, 0, 0,
									0, 0, 0,
									50, 0, 0,
									-50, 50, 0,
									0, 50, 0,
									50, 50, 0,
									-50, -50, 50,
									0, -50, 50,
									50, -50, 50,
									-50, 0, 50,
									0, 0, 50,
									50, 0, 50,
									-50, 50, 50,
									0, 50, 50,
									50, 50, 50 };

								dd2[0] = location[i][0];
								dd2[1] = location[i][1];
								dd2[2] = location[i][2];
								printf("dd2:%f %f %f\n", dd2[0], dd2[1], dd2[2]);
							}
							//////////ループ開始/////////////////////////////////////////////
		//					for (LOOP = 0; LOOP < 2000000; LOOP++) {
							for (i_p = 0; i_p < 3; i_p++)
								pos[0][i_p] = dd2[i_p];
							for (i_a = 0; i_a < 2; i_a++)
								angles[0][i_a] = dd2[3 + i_a];
							for (p = 0; p < update_max; p++)
								S[p] = DBL_MAX;//二乗予測誤差規範 初期化
	//						for (LOOP = 1; LOOP < update_max; LOOP++) {//LOOP:更新の繰り返し回数
								//				Calc_dd(Zy,i,POS,angles,LOOP-1,Scoil_pole,Qm,&S[LOOP]);
							if (Calc_dd_loop(data_buf, i, pos, angles, update_max, &update_num[b], Scoil_pole, Qm, S) != 0) {
								//					printf("ddでか杉でbreak!更新%d回目\n",LOOP);
								//					S[LOOP]=DBL_MAX;
								break;
							}
							//				printf("予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[LOOP]);
//							if(S[LOOP] >= S[LOOP-1])
//								break;
//							I2 = min_index(S, 0, update_max);
//							S1[0] = S[I2];//S1はIAlenだけ確保してある
//						}
							I2 = min_index(S, 0, update_max);
							S1[0] = S[I2];//S1はIAlenだけ確保してある
							//printf("S1:%lf %lf %lf %d IA,k\n",S1[0],Initial_Angle[k][0],Initial_Angle[k][1],k);
							for (i_p = 0; i_p < 3; i_p++)
								dd_stock[b][i_p] = dd2[i_p] = pos[I2][i_p];
							for (i_a = 0; i_a < 2; i_a++)
								dd_stock[b][3 + i_a] = dd2[3 + i_a] = angles[I2][i_a];
							sum2[b] = S1[0];
							sum3[b] = S1[0];

							//						for (j = 0; j < 6; j++)
							//							for (k = 0; k < 1; k++) {
							//								//dd[0]=0;dd[1]=0;dd[2]=0;
							//								//printf("DD:%lf\t%lf\t%lf\t%lf\t%lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
							//								//printf("DD2:%lf\t%lf\t%lf\t%lf\t%lf\n",dd2[0],dd2[1],dd2[2],dd2[3],dd2[4]);
							//								Dfuncg_dd(dummy[i], 1.0, Scoil[j], Z[j], k, u, dd2);
							//								//printf("Z[%d]:%lf\n",k,Z[j][k]);
							//								SDM_2(dummy, dd2, k, Scoil[j], A2[j][k]);
							//								//printf("Z[%d]:%lf\n",k,Z[i][k]);
							//								//printf("A %lf %lf %lf %lf %lf\n",A2[i][k][0],A2[i][k][1],A2[i][k][2],A2[i][k][3],A2[i][k][4]);
							//							}
							//						for (j = 0; j < 5; j++)
							//							d_f[j] = 0.0;
							//						for (u = 0; u < 6; u++)
							//							for (k = 0; k < 1; k++)
							//								for (j = 0; j < 5; j++) {
							//									d_f[j] += -2.0 * (data_buf[i][u] - Qm[u] * Z[u][k]) * (Qm[u] * A2[u][k][j]);
							//									//printf("d_f[%f]",d_f[j]);
							//									//printf("d_buf[%f]",data_buf[k][i]);
							//								}
							//						//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",d_f[0],d_f[1],d_f[2],d_f[3],d_f[4]);
							//						for (k = 0; k < PATERN; k++) {
							//							for (j = 0; j < 5; j++)
							//								dd_buf[k][j] = dd2[j] - alpha[k] * d_f[j];
							//							//printf("dd_buf:%lf:%lf:%lf:%lf:%lf\n",dd_buf[k][0],dd_buf[k][1],dd_buf[k][2],dd_buf[k][3],dd_buf[k][4]);
							//							sum2[k] = 0.0;
							//							sum3[k] = 0.0;
							//							for (u = 0; u < 6; u++)
							//								for (j = 0; j < 1; j++)	
							//									Dfuncg_dd(dummy[i], 1.0, Scoil[u], Z[u], j, u, dd_buf[k]);
							//									//printf("Z[%d]:%lf\n",j,Z[i][j]);
							//									//printf("dd_buf:%lf:%lf\n",dd_buf[k][3],dd_buf[k][4]);
							//									//printf("[%d][%d]\t%f\n",k,j,Q[i]*Z[i][j]);
							//							for (u = 0; u < 6; u++)
							//								for (j = 0; j < 1; j++) {
							//									//printf("[%d][%d]\t%f\t%f\n",k,j,calib_buf[i][j],Q[i]*Z[i][j]);
							//									sum2[k] += pow((data_buf[i][u] - Qm[u] * Z[u][j]), 2.0);
							//									sum3[k] += pow((data_buf[i][u] - Qm[u] * Z[u][j]), 2.0);
							//									//sum2[k] += pow((calib_buf[i][j+(N*calib_point)]- Q[i]*Z[i][j+(N*calib_point)]),2.0);
							//								}
							//							//sum2[k]=sum2[k]/2.0;
							//							//printf("sum2[%d]\t%f\n",k,sum2[k]);
							//						}
							//						MIN = min_index(sum2, 0, PATERN);
							//						for (k = MIN; k < MIN + 1; k++) {
							//							for (j = 0; j < 5; j++)
							//								dd_buf[k][j] = dd2[j] - alpha[k] * d_f[j];
							//							//printf("dd_buf:%lf:%lf:%lf:%lf:%lf\n",dd_buf[k][0],dd_buf[k][1],dd_buf[k][2],dd_buf[k][3],dd_buf[k][4]);
							//							//sum2[k]=0.0;
							//							for (u = 0; u < 6; u++)
							//								for (j = 0; j < 1; j++)
							//									Dfuncg_dd(dummy[i], 1.0, Scoil[u], Z[u], j, u, dd_buf[k]);
							//									//printf("Z[%d]:%lf\n",j,Z[i][j]);
							//									//printf("dd_buf:%lf:%lf\n",dd_buf[k][3],dd_buf[k][4]);
							//									//printf("[%d][%d]\t%f\n",k,j,Q[i]*Z[i][j]);
							//							//printf("d_f:%lf:%lf:%lf:%lf:%lf:%lf\n",pow((data_buf[i][0]- Qm[0]*Z[0][0]),2.0),pow((data_buf[i][1]- Qm[1]*Z[1][0]),2.0),pow((data_buf[i][2]- Qm[2]*Z[2][0]),2.0),pow((data_buf[i][3]- Qm[3]*Z[3][0]),2.0),pow((data_buf[i][4]- Qm[4]*Z[4][0]),2.0),pow((data_buf[i][5]- Qm[5]*Z[5][0]),2.0));
							//							//pow((data_buf[i][u]- Qm[u]*Z[u][j]),2.0);
							//						}
							//						//MIN=min_index(sum2,0,PATERN);
							//						//if(LOOP==0)
							//						//printf("sum2:%lf\n",sum2[MIN]);
							//						for (j = 0; j < 5; j++)
							//							dd2[j] = dd2[j] - alpha[MIN] * d_f[j];
							//						//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",dd[0],dd[1],dd[2],dd[3],dd[4]);
							//						//printf("sum2:%lf\n",pow(sum2[MIN],0.5));
							//						//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",dd2[0]-dd[0],dd2[1]-dd[1],dd2[2]-dd[2],dd2[3],dd2[4]);
							//						//printf("dd:%lf\n",pow((pow((RLP[i][0]-(dd2[0]-dd[0])),2.0)+pow((RLP[i][1]-(dd2[1]-dd[1])),2.0)+pow((RLP[i][2]-(dd2[2]-dd[2])),2.0)),0.5));
							//						if (LOOP == 0)
							//							sum_buf_10000 = sum2[MIN];
							//						//printf("%lf\n",sum2[MIN]);
							//						//printf("%lf\n",sum_buf_10000);
							//						if (LOOP > 999)
							//							if (LOOP % 1000 == 0) {
							//								if((sum2[MIN] / sum_buf_10000) > 0.99999) {
							//									puts("LOOOPEND\n");
							//									break;
							//								}
							//								//printf("%lf\n",sum2[MIN] / sum_buf_10000);
							//								//printf("d_f:%lf:%lf:%lf:%lf:%lf\n",dd2[0]-dd[0],dd2[1]-dd[1],dd2[2]-dd[2],dd2[3],dd2[4]);
							//								//printf("sum2:%lf\n",pow(sum2[MIN],0.5));
							//								sum_buf_10000 = sum2[MIN];
							//							}
							//					}
												/////////////ループ終了////////////////////////////////////////////
							if (zerocorrection == 1)
								printf("d_f[%d]:%lf:%lf:%lf:%lf:%lf\n", i, dd2[0] - dd[0], dd2[1] - dd[1], dd2[2] - dd[2], dd2[3], dd2[4]);
							else
								printf("d_f[%d]:%lf:%lf:%lf:%lf:%lf\n", i, dd2[0], dd2[1], dd2[2], dd2[3], dd2[4]);

							printf("sum2[%d][%d][%d][%d]:%lf\n", i_data_no, i_ch, b, update_num[b] - 1, pow(sum3[b] / NUM_SENT_COIL, 0.5));
							/*
							for(p=0;p<update_max;p++){
							S[p]=DBL_MAX;//予測誤差規範 初期化
							}
							for(j=1;j<update_max;j++){//j:更新の繰り返し回数
							Calc_dd(data_buf,i,pos,angles,j-1,Scoil_pole,Qm,&S[j]);
							//printf("二乗予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[j]);
							if( S[j] > S[j-1]){
							//printf("%d:更新%d回でbreak\n",i,j);
							//printf("二乗予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[j]);
							break;
							}
							}
							I=min_index(S,0,update_max);//更新を繰り返す中で一番誤差が小さかった更新回数 I
							//signal_error[i]=pow(S[I],0.5);//信号誤差 1 x Data_Point_Number
							*/
							//signal_error[i]=pow(sum2[MIN],0.5);
							signal_error[i] = pow(sum3[b] / NUM_SENT_COIL, 0.5);
							// 			printf("%d 予測誤差規範の値::%lf\n",i,S[I]);
							//printf("\n\t%d:%f\n",i,signal_error[i]);
							//printf("%f:%f:%f\n",pos[I][0]-dd[0],pos[I][1]-dd[1],pos[I][2]-dd[2]);
							//		puts("ok1");
							if (U_N == 0) {
								s0 = 1000.0;
								for (j = 0; j < 5; j++)
									//POSAmS[i][j] = pos[I][j];
									POSAmS[i][j] = dd2[j];
							}
							if ((s0 > signal_error[i]) &&  //推定結果が測定領域内にあるか
								(dd2[0] < 70.0 + dd[0]) && (dd2[0] > -70.0 + dd[0]) &&
								(dd2[1] < 70.0 + dd[1]) && (dd2[1] > -70.0 + dd[1]) &&
								(dd2[2] < 70.0 + dd[2]) && (dd2[2] > -70.0 + dd[2])
								//&&((((angles[I][0]*180.0/M_PI) < 40.0) && ((angles[I][0]*180.0/M_PI) > -40.0)) || 
								//((((angles[I][0]*180.0/M_PI) > 140.0) && ((angles[I][0]*180.0/M_PI) < 180.0)) ||
								//(((angles[I][0]*180.0/M_PI) < -140.0) && ((angles[I][0]*180.0/M_PI) > -180.0)))) &&		
								//((((angles[I][1]*180.0/M_PI) < 40.0) && ((angles[I][1]*180.0/M_PI) > -40.0)) || 
								//((((angles[I][1]*180.0/M_PI) > 140.0) && ((angles[I][1]*180.0/M_PI) < 180.0)) ||
								//(((angles[I][1]*180.0/M_PI) < -140.0) && ((angles[I][1]*180.0/M_PI) > -180.0))))
								) {
								//puts("ok");
								//		printf("%f\t%f\t%f\t%f\t%f\t%f\n",signal_error[i],pos[I][0],pos[I][1],pos[I][2],angles[I][0]*180.0/M_PI,angles[I][1]*180.0/M_PI);
								//s0=signal_error[i];
								s0 = pow(sum2[b], 0.5);
								for (j = 0; j < 5; j++)
									//POSAmS[i][j] = pos[I][j];
									POSAmS[i][j] = dd2[j];
								//printf("%d:%d:%f\n",i,U_N,signal_error[i]);	
							}
							//printf("%d:%d:%f\n",i,U_N,s0);
							//printf("%d:%f:%f:%f:%f:%f\n",i,POSAmS[i][0],POSAmS[i][1],POSAmS[i][2],POSAmS[i][3]*180.0/M_PI,POSAmS[i][4]*180.0/M_PI);
		//				}
						//signal_error[i]=pow(sum2[MIN],0.5);
							signal_error[i] = pow(sum3[b] / NUM_SENT_COIL, 0.5);
							//signal_error[i]=s0;
							//printf("%d:%f\n",i,signal_error[i]);
							//printf("%d 予測誤差規範の値::%lf\n",i,signal_error[i]);
			//				if (signal_error[i] > 4.0/*1000000000.0*/) {  ////  L_error_indexの評価はここだよ
							if (signal_error[i] > PARMIT_ERROR_LEVEL /*1000000000.0*/) {  ////  L_error_indexの評価はここだよ
								if (b >= PATERN - 1) { // １つ目のデータ(b==0)か、もしくは前回が40以上で繰り返し、もしくIA初期値にしても大きかったとき
									L_error_index[a] = i; // L_error_indexを記録し
									a++;                // aをインクリメントして
									b = 0;				// 40以上だけど次のデータに行く		
								}
								else {
									b++;
									i--;
								}		//前回が1のとき,もう一回計算
							}
							else {//signal_error[i]が小さかったとき
							 //printf("%d:%f\n",i,signal_error[i]);
							 //printf("ok\n");
								MIN = b;
								sum += signal_error[i];
								//sum +=1.0;
								old_b = b = 1;//b=0でいつでもIAからになるよ 普段は1 Calib_Cdoubleのほうでは0と同じ
								//printf("%lf\n",sum);
							}
							if (b == 0) {
								MIN = min_index(sum3, old_b, PATERN);
								for (i_p = 0; i_p < 3; i_p++)
									dd2[i_p] = dd_stock[MIN][i_p];
								for (i_a = 0; i_a < 2; i_a++)
									dd2[3 + i_a] = dd_stock[MIN][3 + i_a];

								if (zerocorrection == 1)
									printf("d_f[%d]:%lf:%lf:%lf:%lf:%lf\n", i, dd2[0] - dd[0], dd2[1] - dd[1], dd2[2] - dd[2], dd2[3], dd2[4]);
								else
									printf("d_f[%d]:%lf:%lf:%lf:%lf:%lf\n", i, dd2[0], dd2[1], dd2[2], dd2[3], dd2[4]);
								printf("sum2[%d][%d][%d][%d]:%lf\n", i_data_no, i_ch, MIN, update_num[MIN] - 1, pow(sum3[MIN] / NUM_SENT_COIL, 0.5));
								/*
				for(p=0;p<update_max;p++){
				S[p]=DBL_MAX;//予測誤差規範 初期化
				}
				for(j=1;j<update_max;j++){//j:更新の繰り返し回数
				Calc_dd(data_buf,i,pos,angles,j-1,Scoil_pole,Qm,&S[j]);
				//printf("二乗予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[j]);
				if( S[j] > S[j-1]){
				//printf("%d:更新%d回でbreak\n",i,j);
				//printf("二乗予測誤差規範の値,Calc_ddから戻った後メイン:%lf\n",S[j]);
				break;
				}
				}
				I=min_index(S,0,update_max);//更新を繰り返す中で一番誤差が小さかった更新回数 I
				//signal_error[i]=pow(S[I],0.5);//信号誤差 1 x Data_Point_Number
				*/
				//signal_error[i]=pow(sum2[MIN],0.5);
//					signal_error[i] = pow(sum3[MIN] / 6.0, 0.5);
					// 			printf("%d 予測誤差規範の値::%lf\n",i,S[I]);
					//printf("\n\t%d:%f\n",i,signal_error[i]);
					//printf("%f:%f:%f\n",pos[I][0]-dd[0],pos[I][1]-dd[1],pos[I][2]-dd[2]);
					//		puts("ok1");
								if (U_N == 0) {
									s0 = 1000.0;
									for (j = 0; j < 5; j++)
										//POSAmS[i][j] = pos[I][j];
										POSAmS[i][j] = dd2[j];
								}
								if ((s0 > signal_error[i]) &&  //推定結果が測定領域内にあるか
									(dd2[0] < 70.0 + dd[0]) && (dd2[0] > -70.0 + dd[0]) &&
									(dd2[1] < 70.0 + dd[1]) && (dd2[1] > -70.0 + dd[1]) &&
									(dd2[2] < 70.0 + dd[2]) && (dd2[2] > -70.0 + dd[2])
									//&&((((angles[I][0]*180.0/M_PI) < 40.0) && ((angles[I][0]*180.0/M_PI) > -40.0)) || 
									//((((angles[I][0]*180.0/M_PI) > 140.0) && ((angles[I][0]*180.0/M_PI) < 180.0)) ||
									//(((angles[I][0]*180.0/M_PI) < -140.0) && ((angles[I][0]*180.0/M_PI) > -180.0)))) &&		
									//((((angles[I][1]*180.0/M_PI) < 40.0) && ((angles[I][1]*180.0/M_PI) > -40.0)) || 
									//((((angles[I][1]*180.0/M_PI) > 140.0) && ((angles[I][1]*180.0/M_PI) < 180.0)) ||
									//(((angles[I][1]*180.0/M_PI) < -140.0) && ((angles[I][1]*180.0/M_PI) > -180.0))))
									) {
									//puts("ok");
									//		printf("%f\t%f\t%f\t%f\t%f\t%f\n",signal_error[i],pos[I][0],pos[I][1],pos[I][2],angles[I][0]*180.0/M_PI,angles[I][1]*180.0/M_PI);
									//s0=signal_error[i];
									s0 = pow(sum2[MIN], 0.5);
									for (j = 0; j < 5; j++)
										//POSAmS[i][j] = pos[I][j];
										POSAmS[i][j] = dd2[j];
									//printf("%d:%d:%f\n",i,U_N,signal_error[i]);	
								}
								signal_error[i] = pow(sum3[MIN] / NUM_SENT_COIL, 0.5);
								//printf("%d:%d:%f\n",i,U_N,s0);
								//printf("%d:%f:%f:%f:%f:%f\n",i,POSAmS[i][0],POSAmS[i][1],POSAmS[i][2],POSAmS[i][3]*180.0/M_PI,POSAmS[i][4]*180.0/M_PI);
								old_b = b = 1;
							}
						}
					}
					////////////////////////////////////////////////////////
					//              !!!! Lines loop END !!!!!             //
					////////////////////////////////////////////////////////
					//for(i=0;i<Lines;i++)
					//{printf("\n%d:%f\n",i,signal_error[i]);
					//}
					////////////////////-----------角度誤差計算と-pi～piへの丸め
					/*
					for(i=0;i<Lines;i++){
					angle_wrap2(&POSAmS[i][3]);angle_wrap2(&POSAmS[i][4]);
					if(fabs((double)(POSAmS[i][3]+fabs(POSAmS[i][4]) > M_PI*7/6))){
					POSAmS[i][3]=-POSAmS[i][3]+M_PI;POSAmS[i][4]+=M_PI;
					}
					angle_wrap2(&POSAmS[i][3]);angle_wrap2(&POSAmS[i][4]);
					*/
					//////////////////---------------------------------------------------
					for (i = 0; i < Lines; i++)
						if (signal_error[i] < 40)
							printf("POSAmS:%d:\t%lf\t%lf\t%lf\t%lf\t%lf\t%+e \n", i, POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
						else
							printf("POSAmS:%d:\t%lf\t%lf\t%lf\t%lf\t%lf\t%+e \t failed\n", i, POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
					mean_signal_error = sum / (Lines - a);
					printf("a:%d\n", a);
					printf("mean_signal_error:%lf  failed %d points\n", mean_signal_error, a);
					//for(i=0;i<Lines;i++)
					//{printf("\n%d:%f\n",i,signal_error[i]);
					//}
					//scanf("%d",&i);
					//////____位置誤差の計算______evaluation==1のときのみ____________________________________
					for (i = 0; i < Lines; i++) {
						printf("RLP[%d]:%lf %lf %lf \n", i, RLP[i][0], RLP[i][1], RLP[i][2]);
						printf("POSAmS[%d]:%lf %lf %lf \n", i, POSAmS[i][0], POSAmS[i][1], POSAmS[i][2]);
					}
					POSerrormean = 0.0;
					a = 0;
					for (j = 0; j < 3; j++)
						res[j] = 0.0;
					for (i = 0; i < Lines; i++) {
						POSerror[i] = 0.0;
						for (j = 0; j < 3; j++) {
							printf("%f", RLP[i][j]);
							POSerror[i] += pow((POSAmS[i][j] - RLP[i][j]), 2.0);  //理論上の受信コイル位置情報からの、計算した受信コイルの位置の	
						}
						POSerror[i] = pow(POSerror[i], 0.5);
						if (L_error_index[a] == i) {
							printf("POSerror[%d]:%lf[mm]\tBut estimation failed at this point!\n", i, POSerror[i]);
							a++;
						}
						else {
							POSerrormean += POSerror[i];
							printf("POSerror[%d]:%lf[mm]\n", i, POSerror[i]);
							for (j = 0; j < 3; j++)
								res[j] += (POSAmS[i][j] - RLP[i][j]);
						}
					}
					for (j = 0; j < 3; j++) {
						//res[j] /= (Lines-a);
						res[j] = dd[j];
						printf("res_%lf\n", res[j]);
					}
					//printf("anyway...enter some number...");
					//scanf("%d",&ss);
					res_error_mean = 0.0;
					a = 0;
					for (i = 0; i < Lines; i++) {
						res_error[i] = 0.0;
						for (j = 0; j < 3; j++)
							res_error[i] += pow((POSAmS[i][j] - RLP[i][j] - res[j]), 2.0);
						res_error[i] = pow(res_error[i], 0.5);
						if (L_error_index[a] == i) {
							printf("POSerror[%d]:%lf[mm]\tBut estimation failed at this point!\n", i, res_error[i]);
							a++;
							max_buf[i] = 0.0;
						}
						else {
							res_error_mean += res_error[i];
							printf("res_error[%d]:%lf[mm]\n", i, res_error[i]);
							max_buf[i] = res_error[i];
							angle_av[0] += POSAmS[i][3];
							angle_av[1] += POSAmS[i][4];
						}
					}
					// POSerrormean = sum(POSerror)/numel(POSerror);
					POSerrormean /= (Lines - a);
					res_error_mean /= (Lines - a);
					angle_av[0] /= (Lines - a);
					angle_av[1] /= (Lines - a);
					printf("POSerrormean:%lf[mm] 点数:%d\n", POSerrormean, Lines - a);
					printf("res_errormean:%lf[mm] 点数:%d\n", res_error_mean, Lines - a);
					MAX = max_index(max_buf, 0, Lines);
					//////____位置誤差の計算______evaluation==1のときのみ____ここまで____________________________
				}
				///////////////    データの書き出し---------------------
				////位置のデータをプロット用に書き出す
				sprintf(out_fldir, "..\\..\\data\\%s\\%s\\S%s_D%s", DATA_DATE, Data_dir, Sc_dir, Data_dir);
				_mkdir(out_fldir);
				sprintf(out_flname, "..\\..\\data\\%s\\%s\\S%s_D%s\\S%s_D%s_%d_ch%d_POSest_result.data", DATA_DATE, Data_dir, Sc_dir, Data_dir, Sc_dir, Data_dir, file_no, i_ch);
				printf("%s\n", out_flname);
				if ((fp = fopen(out_flname, "w+")) == NULL) {
					printf("file open error @file write POSest_result.data\n");
					return EXIT_FAILURE;
				}
				a = 0;
				for (i = 0; i < Lines; i++) {
					printf("Lines = %d\n", i);
					if (zerocorrection == 1)
						if (L_error_index[a] == i) {
							fprintf(fp, "# %lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0] - dd[0], POSAmS[i][1] - dd[1], POSAmS[i][2] - dd[2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
							a++;
						}
						else
							fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0] - dd[0], POSAmS[i][1] - dd[1], POSAmS[i][2] - dd[2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
					else
						if (L_error_index[a] == i) {
							fprintf(fp, "# %lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
							a++;
						}
						else
							fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
				}
				fclose(fp);
				printf("good\n");
				//////////////////////////////////
				//////////------evaluation==1のとき//////////////////////////////////////////////////////////
				if (evaluation == 1) {
					//////---POSerror,ANGLEerror,Signal_errorを書き出す
					//sprintf(out_flname,"./evaluation_results/calibed_by_%s",ScoilDirname);
					//mkdir(out_flname,0777);
					//sprintf(out_flname,"./evaluation_results/calibed_by_%s/%s",ScoilDirname,Scoilfilename);
					//mkdir(out_flname,0777);
					//sprintf(out_fldir, "..\\..\\data\\%s\\S%s_D%s", Data_dir, Sc_dir, Data_dir);
					//mkdir(out_fldir, 0777);
					sprintf(out_flname, "..\\..\\data\\%s\\%s\\S%s_D%s\\S%s_D%s_ch%d_POS_angle.data", DATA_DATE, Data_dir, Sc_dir, Data_dir, Sc_dir, Data_dir, i_ch);
					printf("%s\n", out_flname);
					if ((fp = fopen(out_flname, "w+")) == NULL) {
						printf("file open error @file write POSest_result.data\n");
						return EXIT_FAILURE;
					}
					a = 0;
					for (i = 0; i < Lines; i++) {
						if (zerocorrection == 1)
							if (L_error_index[a] == i) {
								fprintf(fp, " %lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0] - dd[0], POSAmS[i][1] - dd[1], POSAmS[i][2] - dd[2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
								a++;
							}
							else
								fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0] - dd[0], POSAmS[i][1] - dd[1], POSAmS[i][2] - dd[2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
						else
							if (L_error_index[a] == i) {
								fprintf(fp, " %lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
								a++;
							}
							else
								fprintf(fp, "%lf\t%lf\t%lf\t%lf\t%lf\t%+e\n", POSAmS[i][0], POSAmS[i][1], POSAmS[i][2], POSAmS[i][3], POSAmS[i][4], signal_error[i]);
					}
					fclose(fp);
					sprintf(out_flname, "..\\..\\data\\%s\\%s\\S%s_D%s\\S%s_D%s_ch%d_pos_angle_signalerror.data", DATA_DATE, Data_dir, Sc_dir, Data_dir, Sc_dir, Data_dir, i_ch);
					printf("%s\n", out_flname);
					if ((fp = fopen(out_flname, "w+")) == NULL) {
						printf("file open error @file write pos_angle_signalerror.data\n");
						return EXIT_FAILURE;
					}
					a = 0;
					fprintf(fp, "POS_e\tres_e\tsignal_e\n");
					for (i = 0; i < Lines; i++) {
						if (L_error_index[a] == i) {
							fprintf(fp, "%lf\t%lf\t%+e\n", POSerror[i], res_error[i], signal_error[i]);
							a++;
						}
						else fprintf(fp, "%lf\t%lf\t%+e\n", POSerror[i], res_error[i], signal_error[i]);
					}
					fclose(fp);
					/////  POSerrormean,mean_Signal_error
					sprintf(out_flname, "..\\..\\data\\%s\\%s\\S%s_D%s\\S%s_D%s_ch%d_error_mean.data", DATA_DATE, Data_dir, Sc_dir, Data_dir, Sc_dir, Data_dir, i_ch);
					if ((fp = fopen(out_flname, "w+")) == NULL) {
						printf("file open error @file write error_mean.txt\n");
						return EXIT_FAILURE;
					}
					fprintf(fp, "POS_e\tres_e\tsignal_e\n");
					fprintf(fp, "%lf\t%lf\t%lf\n", POSerrormean, res_error_mean, mean_signal_error);
					fprintf(fp, "\nmean value excluded data at %d large error point\n", a);
					//// signal errorを中心データで基準化
					if (DP == 5)
						fprintf(fp, "signal error regulared by signal:%lf\n", mean_signal_error / pow(pow(data_buf[62][0], 2.0) + pow(data_buf[62][1], 2.0) + pow(data_buf[62][2], 2.0) + pow(data_buf[62][3], 2.0) + pow(data_buf[62][4], 2.0) + pow(data_buf[62][5], 2.0) + pow(data_buf[62][6], 2.0) + pow(data_buf[62][7], 2.0), 0.5));
					if (DP == 3)
						fprintf(fp, "signal error regulared by signal:%lf\n", mean_signal_error / pow(pow(data_buf[14][0], 2.0) + pow(data_buf[14][1], 2.0) + pow(data_buf[14][2], 2.0) + pow(data_buf[14][3], 2.0) + pow(data_buf[14][4], 2.0) + pow(data_buf[14][5], 2.0) + pow(data_buf[14][6], 2.0) + pow(data_buf[14][7], 2.0), 0.5));
					if (DP == 7)
						fprintf(fp, "signal error regulared by signal:%lf\n", mean_signal_error / pow(pow(data_buf[3][0], 2.0) + pow(data_buf[3][1], 2.0) + pow(data_buf[3][2], 2.0) + pow(data_buf[3][3], 2.0) + pow(data_buf[3][4], 2.0) + pow(data_buf[3][5], 2.0) + pow(data_buf[3][6], 2.0) + pow(data_buf[3][7], 2.0), 0.5));
					if (DP == 9)
						fprintf(fp, "signal error regulared by signal:%lf\n", mean_signal_error / pow(pow(data_buf[4][0], 2.0) + pow(data_buf[4][1], 2.0) + pow(data_buf[4][2], 2.0) + pow(data_buf[4][3], 2.0) + pow(data_buf[4][4], 2.0) + pow(data_buf[4][5], 2.0) + pow(data_buf[4][6], 2.0) + pow(data_buf[4][7], 2.0), 0.5));
					if (DP == 75)
						fprintf(fp, "signal error regulared by signal:%lf\n", mean_signal_error / pow(pow(data_buf[37][0], 2.0) + pow(data_buf[37][1], 2.0) + pow(data_buf[37][2], 2.0) + pow(data_buf[37][3], 2.0) + pow(data_buf[37][4], 2.0) + pow(data_buf[37][5], 2.0) + pow(data_buf[37][6], 2.0) + pow(data_buf[37][7], 2.0), 0.5));
					fprintf(fp, "%lf\t%lf\t%lf\n", res[0], res[1], res[2]);
					fclose(fp);
				}
				//////////------evaluation==1のとき--------------ここまで------------
				printf("good2\n");
				if (evaluation == 1)
					fprintf(rslt_fp, "%s_%s_%d_%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Sc_dir, Data_dir, file_no, i_ch, POSerrormean, res_error_mean, mean_signal_error, res[0], res[1], res[2], res_error[MAX], signal_error[MAX], angle_av[0], angle_av[1]);
				else
					printf("MAX=%d\n", MAX);
				fprintf(rslt_fp, "%s_%s_%d_%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", Sc_dir, Data_dir, file_no, i_ch, mean_signal_error, res[0], res[1], res[2], res_error[MAX], signal_error[MAX], angle_av[0], angle_av[1]);
				printf("good3\n");
				free_matrix(data_buf, Lines + 1);
				free_matrix(pos, update_max);
				free_matrix(angles, update_max);
				free(signal_error2);
				free(mean_signal_error2);
				free(signal_error);
				free_matrix(POSAmS, Lines);
				free(S);
				free(S1);
				free(L_error_index);
				free_matrix(Scoil_pole, NUM_SENT_COIL);
				free_matrix(Initial_Angle, IAlen);
				if (evaluation == 1) {
					free(POSerror);
					free(res_error);
					free(max_buf);
					free(dd);
					free(dd2);
					free_matrix(RLP, Lines);
					free_matrix(dummy, Lines);
					free_matrix(dd_stock, PATERN);
					free_matrix(Z, NUM_SENT_COIL);
				}
				//		fclose(fp);
			}
			free_matrix(U_B, 9 * 25);
			free_matrix(Scoil, NUM_SENT_COIL);
			free_matrix(dd_read, FILENUM);
			free(Qm);
			fclose(rslt_fp);
		}
	}

	return 0;
}
