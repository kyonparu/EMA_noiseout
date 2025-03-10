// sumi2.cpp : コンソール アプリケーションのエントリ ポイントを定義します。
//

#include "stdafx.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<float.h>
/*
#define _USE_MATH_DEFINES
*/

#define M_PI 3.141592653589793
//#define CL 60.0			//コイル長[mm]
#define CL 58.0
#define L  337.5		//座標中心からコイルの中心までの長さ[mm]
#define PATERN 8		//収束因子の数
#define FILENUM 1
// #define FILENUM 9		//校正用データセットの数
//#define calib_point 27	//データセットの校正点の数
#define calib_point 189
#define	GL	100			//格子の一辺の長さ[mm]
#define Switch 0		//受信コイルの位置を校正する場合は0
#define r_Switch 0		//オフセットを計算する場合は0
#define w_Switch 1		//ゲインに基づく重みづけをする場合は１
#define rr_Switch 1		//オフセットをデータセットごとに計算する場合は0
#define START_CH_NO 1	//受信コイル最小No.
#define END_CH_NO 12	//受信コイル最大No.
#define NUM_RES_COIL (END_CH_NO - START_CH_NO + 1) //受信コイル数
#define REF_CH_NO 2
#define NUM_SEN_COIL 8 //送信コイル数


double **new_matrix(int m, int n)
{
  int i,j;
  double **a;
  a = (double **)malloc(sizeof(double *) * m);
  for (i = 0; i < m; i++)
	  a[i] = (double *)malloc(sizeof(double) * n);

for(i=0;i<m;i++){
	for(j=0;j<n;j++){
		a[i][j]=0;
		}
	}
  return a;
}

int min_index(double *S,int i1,int i2){
int i,j;
double C;
j=i1;
C=S[i1];
for(i=i1+1;i<i2;i++){
	if(S[i] < C){
//printf("%d\n",i);
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

void Dfuncg_dd(double **X,double qm,double *s,double *G,int T2,int I,double *dd,int N){
int i,j,k;
double Q,Rp[3],Rn[3];
double lenRp=0,lenRn=0,Yp[3],Yn[3],Y[3],CC1,CS1,Si1,CC2,CS2,Si2,CCX,CSX,SiX,CCY,CSY,SiY,CCZ,CSZ,SiZ,CC,CS,Si;
lenRp=0.0;
lenRn=0.0;
//U0 = 1.26/(10^6);      //空気の透磁率≒真空の透磁率＝1.26×10^-6

CC = cos(dd[3]/180.0*M_PI)*cos(dd[4]/180.0*M_PI);
CS = cos(dd[3]/180.0*M_PI)*sin(dd[4]/180.0*M_PI);
Si = -sin(dd[3]/180.0*M_PI);

CC2=CC;
CS2=CS;
Si2=Si;



    Q=qm/(4.0*M_PI*(1.26/pow(10.0,4.0)));


//S:1x6 coilNo番目の両極の座標

    Rp[0] = X[T2][0]+dd[0] - s[0];             //p極から点X(今の受信信号の位置)までのベクトル
    Rp[1] = X[T2][1]+dd[1] - s[1];
    Rp[2] = X[T2][2]+dd[2] - s[2];
    Rn[0] = X[T2][0]+dd[0] - s[3];             
    Rn[1] = X[T2][1]+dd[1] - s[4];
    Rn[2] = X[T2][2]+dd[2] - s[5];            //n極から点Xまでのベクトル


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
//    printf("Yn[%d]:%lf ",i,Yn[i]);
//    printf("Y[%d]:%lf\n",i,Y[i]);
    }
////-------------

G[T2+N] = Y[0]*CC2 + Y[1]*CS2 + Y[2]*Si2;

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

void SDM_1(double *Scoil_pole,double *dd,double **POS,double *A,int a)
{

int i,j,k;
double CC,CS,Si,Rp[3],Rn[3],lenRp,lenRn,roRp3x1,roRp3x2,roRp3x3,roRn3x1,roRn3x2,roRn3x3;
double g;

g=1.0/(4.0*M_PI*(1.26/pow(10.0,4.0)));

CC = cos(dd[3]/180.0*M_PI)*cos(dd[4]/180.0*M_PI);
CS = cos(dd[3]/180.0*M_PI)*sin(dd[4]/180.0*M_PI);
Si = sin(dd[3]/180.0*M_PI);
 //方向を表す単位ベクトルの成分 Siだけ符号が逆にしてある


//for(i=0; i < NUM_SEN_COIL; i++){    
    for(j=0;j<3;j++){
    Rp[j]=0.0;
    Rn[j]=0.0;
    Rp[j] = POS[a][j]+dd[j] - Scoil_pole[j];//i番目送信コイル正極から受信コイル位置までのベクトル
    Rn[j] = POS[a][j]+dd[j] - Scoil_pole[j+3];//
    }

lenRp=0.0;
lenRn=0.0;

    for(j=0;j<3;j++){
    lenRp=lenRp+pow(Rp[j],2);
    lenRn=lenRn+pow(Rn[j],2);   //ベクトルの長さ
    }

    lenRp=pow(lenRp,0.5);
    lenRn=pow(lenRn,0.5);

	//printf("lenRp %lf\n",lenRp);
	//printf("lenRn %lf\n",lenRn);
	


    roRp3x1 = -3.0 * lenRp * Rp[0];  //偏微分らしいけどなんでこうなるの？やまさんの論文9ページ
    roRp3x2 = -3.0 * lenRp * Rp[1];
    roRp3x3 = -3.0 * lenRp * Rp[2];
    roRn3x1 = -3.0 * lenRn * Rn[0];
    roRn3x2 = -3.0 * lenRn * Rn[1];
    roRn3x3 = -3.0 * lenRn * Rn[2];

	/*
	printf("roRp3x1 %lf\n",roRp3x1);
	printf("roRp3x2 %lf\n",roRp3x2);
	printf("roRp3x3 %lf\n",roRp3x3);
	printf("roRn3x1 %lf\n",roRn3x1);
	printf("roRn3x2 %lf\n",roRn3x2);
	printf("roRn3x3 %lf\n",roRn3x3);
	*/



	A[0] = g * (( -1.0*( pow(lenRp,3) + Rp[0]*roRp3x1)/ pow(lenRp,6))*CC
			+( -1.0*(Rp[1]*roRp3x1)/pow(lenRp,6))*CS
				+((Rp[2]*roRp3x1)/pow(lenRp,6))*Si);

	A[1] = g * (( -1.0*(Rp[0]*roRp3x2)/pow(lenRp,6))*CC
			+ ( -1.0*( pow(lenRp,3) + Rp[1]*roRp3x2)/ pow(lenRp,6))*CS
				+((Rp[2]*roRp3x2)/pow(lenRp,6))*Si);

	A[2] = g * (( -1.0*(Rp[0]*roRp3x3)/pow(lenRp,6))*CC
			+(-1.0*(Rp[1]*roRp3x3)/pow(lenRp,6))*CS
				+ (( pow(lenRp,3) + Rp[2]*roRp3x3)/ pow(lenRp,6))*Si);
			
	A[3] = g * (( ( pow(lenRn,3) + Rn[0]*roRn3x1)/ pow(lenRn,6))*CC
			+( (Rn[1]*roRn3x1)/pow(lenRn,6))*CS
				+(-1.0*(Rn[2]*roRn3x1)/pow(lenRn,6))*Si);

	A[4] = g * (( (Rn[0]*roRn3x2)/pow(lenRn,6))*CC
			+ (( pow(lenRn,3) + Rn[1]*roRn3x2)/ pow(lenRn,6))*CS
				+(-1.0*(Rn[2]*roRn3x2)/pow(lenRn,6))*Si);

	A[5] = g * (( (Rn[0]*roRn3x3)/pow(lenRn,6))*CC
			+((Rn[1]*roRn3x3)/pow(lenRn,6))*CS
				+ (-1.0*( pow(lenRn,3) + Rp[2]*roRn3x3)/ pow(lenRn,6))*Si);

  


//}





}

void SDM_2(double **POS,double *dd,int a,double *Scoil_pole,double *A){

int i,j,k,I;
double Q;
double estG0[6];
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


}

void angle_wrap2(double *a){
	while(*a <= -M_PI+0.6)
	*a = *a + 2.0*M_PI;
	
	while(*a > M_PI+0.6)
		*a = *a - 2.0*M_PI;
}/////表示するときはこっちね

void weight(double **buf,double **w)
{
	int i,N,max_n,min_n;
	double d_buf[FILENUM][NUM_SEN_COIL],d_sum;

	for(i=0; i < NUM_SEN_COIL; i++)
	{
		d_sum=0.0;

		for(N=0;N<FILENUM;N++)
		{
			min_n=min_index(buf[i],N*calib_point,(N+1)*calib_point);
			max_n=max_index(buf[i],N*calib_point,(N+1)*calib_point);
			d_buf[N][i]=buf[i][max_n]-buf[i][min_n];
			printf("d:%d:%d:%lf\n",max_n,min_n,d_buf[N][i]);
			d_sum += d_buf[N][i];
		}
		
		for(N=0;N<FILENUM;N++)
		{
			w[N][i]= d_buf[N][i] / d_sum;
			printf("w:%lf\n",w[N][i]);
		}
	}




}

int main(void)
{

	double **calib_buf;
	FILE *in_fp[FILENUM],*fp_res,*fp_Sp,*fp_GP,*fp_Sc,*fp_Sc_in;
	

	double **tens;
	double **NOI;
	double **Z;
	double **Z2;



	int cal_ch_no, i, i_ch, i_fnum, i_send, j,k,a,s;
	//,N;
	int PP,PP2,PP3;
	int count1=0,count2=0;
	double sum_sum3;
	double e[NUM_SEN_COIL][3];//変換行列
	
	/***** hlm系 *****/
	//double SL[6][3]={  0.0,0.0,-1.0  ,  0.0,-1.0,0.0  ,  -1.0,0.0,0.0  ,  0.0,0.0,1.0   ,  0.0,1.0,0.0  ,  1.0,0.0,0.0 };
	//コイル角度：初期位置
	//double Scoil_angles[6][2] = {0.0       ,0.0        ,0.0       ,45.0     ,0.0        ,-90.0      ,45.0       ,-90.0    ,-90.0    ,0.0        ,135.0      ,0.0};
	//double Scoil_angles[6][2]={-85.000000,-165.000000,180.000000,45.000000,-200.000000,-110.000000,-250.000000,25.000000, 20.000000,-305.000000,-100.000000,200.000000};//送信コイル角度
	/***** md4系 *****/
	/*
	double SL[6][3] = { 0.0,-2.0/sqrt(6.0),1.0/sqrt(3.0)  ,
		                1.0/sqrt(2.0),-1.0/sqrt(6.0),-1.0/sqrt(3.0)  ,   
						1.0/sqrt(2.0),1.0/sqrt(6.0),1.0/sqrt(3.0)  ,  
						0.0,2.0/sqrt(6.0),-1.0/sqrt(3.0)   ,  
						-1.0/sqrt(2.0),1.0/sqrt(6.0),1.0/sqrt(3.0)  ,  
						-1.0/sqrt(2.0),-1.0/sqrt(6.0),-1.0/sqrt(3.0) };//md4系
						*/
	double SL[NUM_SEN_COIL][3] = { 0.0, 2.0 * sqrt(6.0) / 9.0, 1.0 / sqrt(3.0),
		0.0, -2.0 * sqrt(6.0) / 9.0, 1.0 / sqrt(3.0),
		2.0 * sqrt(6.0) / 9.0, 0.0, 1.0 / sqrt(3.0),
		-2.0 * sqrt(6.0) / 9.0, 0.0, 1.0 / sqrt(3.0),
		-2.0 * sqrt(3.0) / 9.0, -2.0 * sqrt(3.0) / 9.0, 1.0 / sqrt(3.0),
		2.0 * sqrt(3.0) / 9.0, -2.0 * sqrt(3.0) / 9.0, 1.0 / sqrt(3.0),
	    -2.0 * sqrt(3.0) / 9.0, 2.0 * sqrt(3.0) / 9.0, 1.0 / sqrt(3.0),
	    2.0 * sqrt(3.0) / 9.0, 2.0 * sqrt(3.0) / 9.0, 1.0 / sqrt(3.0)};//md4系
	/* 初期状態 //
	double Scoil_angles[6][2] = { -35.2644 ,   -150.0,
		                            0.0    ,      0.0,
								  -35.2644 ,    150.0,
								    0.0    ,    120.0,
								   35.2644 ,     90.0,
								   54.7356 ,    150.0 };
	*/
	double Scoil_angles[NUM_SEN_COIL][2] = {-90.0, 0.0,
		                          -80.0, 345.0,
								  -80.0, 95.0,
								  -45.0, 240.0,
								  -70.0, 355.0,
		                          -55.0, 135.0,
		                          -15.0, 320.0,
								  -85.0, 355.0};
	
	double **Scoil_pole;
	double **S_p,S_c[NUM_SEN_COIL][6],S_l[NUM_SEN_COIL],S_b[NUM_SEN_COIL][3];
	double dd[FILENUM][5],d_f[NUM_SEN_COIL][6],d_f2[NUM_SEN_COIL],d_f3[FILENUM][5],A[NUM_SEN_COIL][calib_point*FILENUM][6],A2[FILENUM][NUM_SEN_COIL][calib_point*FILENUM][5];
	int index,MIN, ix, iy, iz;
	double sum[NUM_SEN_COIL],sum2[PATERN],sum3[NUM_SEN_COIL][FILENUM][PATERN],sum3_buf[NUM_SEN_COIL][FILENUM],sum4[FILENUM],sum4_buf[FILENUM];
	double Q[NUM_SEN_COIL],q_a[NUM_SEN_COIL][6],q_b[NUM_SEN_COIL][6],q_a2[FILENUM][NUM_SEN_COIL][5],q_b2[FILENUM][NUM_SEN_COIL][5],q_c[NUM_SEN_COIL],q_d[NUM_SEN_COIL],q_dd1[NUM_SEN_COIL][6],q_dd2[FILENUM][NUM_SEN_COIL][5];
	double theta[NUM_SEN_COIL][2];
	double alpha[PATERN];
	double S_p_buf[NUM_SEN_COIL][PATERN][6],dd_buf[PATERN][FILENUM][5];
	double d_sum,sum_buf=0.0,sum_buf_10000 = 1.0;
	double sns[FILENUM][NUM_SEN_COIL];
	char Scoil_fldir[256];
	char Scoil_flname[256];
	char Sc_dir[256];
	char flname[256];
	char d_name[256];
	char d_dir[256];
	char di_dir[FILENUM][256];
	double theta2[NUM_SEN_COIL][2];
	double r[FILENUM][NUM_SEN_COIL],r2[FILENUM][NUM_SEN_COIL];
	double **w;
	int s_flg=1, ch_sel_flg;
	int LOOP_MAX;
	int aiueo, start_ch_no, end_ch_no, ref_ch_no;
	
	//memory_location//

	calib_buf=new_matrix(NUM_SEN_COIL,calib_point*FILENUM);
	tens=new_matrix(calib_point,3);
	Z=new_matrix(NUM_SEN_COIL,calib_point*FILENUM);
	Z2=new_matrix(NUM_SEN_COIL,calib_point);
	Scoil_pole=new_matrix(NUM_SEN_COIL,6);
	S_p=new_matrix(NUM_SEN_COIL,6);
	w=new_matrix(FILENUM,NUM_SEN_COIL);

	//////////////////
	// Assignment data //
	start_ch_no = START_CH_NO;
	end_ch_no = END_CH_NO;
	ref_ch_no = REF_CH_NO;
	/////////////////////

	for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
		printf("put input file dir(No. %d)\n", i_fnum + 1);
		scanf("%s", di_dir[i_fnum]); //入力ファイルディレクトリの入力
	}

	printf("put res-file dir\n");
	scanf("%s", d_dir);
//	printf("%s\n", d_dir);
	printf("input_Sc_file? : y = 1_n = 2\n");//校正済みの送信コイル位置を読み込むかの判断
	scanf("%d", &s_flg);

	if (s_flg == 1) {
		printf("put Scoil file dir\n");
		scanf("%s", Sc_dir);
		printf("put Scoil file channel no.\n");
		scanf("%d", &cal_ch_no);
	}

	printf("enter_receiver_coil's_angle_y=1\n");
	scanf("%d",&aiueo);

	printf("Do you want to select the start channel No. & the end channel No.? : y = 1_n = 2\n");//校正済みの送信コイル位置を読み込むかの判断
	scanf("%d", &ch_sel_flg);

	if (ch_sel_flg == 1) {
		printf("Please input the start channel No.\n");
		scanf("%d", &start_ch_no);
		printf("Please input the end channel No.\n");
		scanf("%d", &end_ch_no);
		printf("Please input the refference channel No.\n");
		scanf("%d", &ref_ch_no);
	}
	
	//pre_setting////////////////////////////////////////////////////////
	for (i = 0; i < NUM_SEN_COIL; i++)
		for (j = 0; j < 3; j++)
			SL[i][j] = SL[i][j] * L;//送信コイルの中心座標
	
	for (i = 0; i < NUM_SEN_COIL; i++) {
		e[i][0] = cos(Scoil_angles[i][0] / 180.0 * M_PI) * cos(Scoil_angles[i][1] / 180.0 * M_PI);
		e[i][1] = cos(Scoil_angles[i][0] / 180.0 * M_PI) * sin(Scoil_angles[i][1] / 180.0 * M_PI);
		e[i][2] = -sin(Scoil_angles[i][0] / 180.0 * M_PI);
		
		for (j = 0; j < 3; j++)
			Scoil_pole[i][j] = SL[i][j] - CL / 2.0 * e[i][j];//s極（頭）の座標
		for (j = 3; j < 6; j++)
			Scoil_pole[i][j] = SL[i][j - 3] + CL / 2.0 * e[i][j - 3]; //n極（頭）の座標
	}
	
	///////校正済みデータの読み込み////////////////////////////////////
	if (s_flg == 1) {
		sprintf(flname, "..\\..\\data\\%s\\c-mean\\out_ch%d_Scoil.data", Sc_dir, cal_ch_no);
		fp_Sc_in = fopen(flname, "r");
		if (fp_Sc_in == NULL) {
			printf("file_open_error_at_Sc_in\n");
			return EXIT_FAILURE;
		}
		for (i = 0; i < NUM_SEN_COIL; i++)
			fscanf(fp_Sc_in, "%lf%lf%lf%lf%lf%lf", &Scoil_pole[i][0], &Scoil_pole[i][1], &Scoil_pole[i][2], &Scoil_pole[i][3], &Scoil_pole[i][4], &Scoil_pole[i][5]);
	}

	for (i_ch = start_ch_no; i_ch <= end_ch_no; i_ch++) {
		if (i_ch == ref_ch_no)
			continue;
		/////file_open_and_input//////////////////////////////////////////
		////校正データセットを読み込む////////////////////////////////////
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
			sprintf(flname, "..\\..\\data\\%s\\c-mean\\out_ch%d.data", di_dir[i_fnum], i_ch); //校正用データセットの名前
//			printf("%s\n", flname);
			in_fp[i_fnum] = fopen(flname, "r");
			if (in_fp[i_fnum] == NULL) {
				printf("file_open_error_at_%d\n", i_fnum);
				return EXIT_FAILURE;
			}
		}
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			for (i = 0; i < calib_point; i++)
				for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
					if (fscanf(in_fp[i_fnum], "%lf", &calib_buf[i_send][i + i_fnum * calib_point]) != 1) {
						printf("file_input_error_at_%d\n", i_fnum);
						break;
					}
	//確認用
	//for(i=0;i<calib_point*FILENUM;i++)
	//{
	//	printf("[%d] %f %f %f %f %f %f\n",i,calib_buf[0][i],calib_buf[1][i],calib_buf[2][i],calib_buf[3][i],calib_buf[4][i],calib_buf[5][i]);
	//}
	////

	///////////////////////////////////////////////////////////////////


	///out_put_file///////////////////////////////////////////////////

	//printf("put_res-file_name\n");//結果のファイルの名前の入力
	//scanf("%s",Scoil_flname);
		sprintf(Scoil_fldir, "..\\..\\data\\%s\\c-mean", d_dir);
		//		printf("Scoil_fldir = %s\n", Scoil_fldir);
		sprintf(Scoil_flname, "..\\..\\data\\%s\\c-mean\\out_ch%d", d_dir, i_ch);
		sprintf(flname, "%s\\res_out_ch%d.data", Scoil_fldir, i_ch);//校正結果ファイル
		//		printf("flname = %s\n", flname);
		fp_res = fopen(flname,"w");
		if (fp_res == NULL) {
			printf("file_open_error_at_fp_res\n");
			return EXIT_FAILURE;
		}
		sprintf(flname, "%s_Sp.data", Scoil_flname);//校正結果ファイル2
		printf("flname = %s\n", flname);
		fp_Sp = fopen(flname,"w");
		if (fp_Sp == NULL) {
			printf("file_open_error_at_fp_Sp\n");
			return EXIT_FAILURE;
		}
		sprintf(flname, "%s_GP.data", Scoil_flname);//校正ゲイン位置ファイル
		fp_GP = fopen(flname,"w");
		if (fp_GP == NULL) {
			printf("file_open_error_at_fp_GP\n");
			return EXIT_FAILURE;
		}
		sprintf(flname, "%s_Scoil.data", Scoil_flname);//校正磁極位置ファイル
		fp_Sc = fopen(flname,"w");
		if (fp_Sc == NULL) {
			printf("file_open_error_at_fp_Sc\n");
			return EXIT_FAILURE;
		}
	
		/////校正点（格子点）の用意///////////////////////////////////////////////
		if (calib_point == 27) {
			PP = 1;
			PP2 = 1;
			PP3 = 1;
		}
		if (calib_point == 125) {
			PP = 2;
			PP2 = 2;
			PP3 = 2;
		}
		if (calib_point == 27 || calib_point == 125) {
			a = 0;
			for (i = -PP3; i <= PP3; i++)
				for (j = -PP2; j <= PP2; j++)
					for(k = -PP; k <= PP; k++) {
						tens[a][0] = ((double)k / PP) * GL / 2;
						tens[a][1] = ((double)j / PP2) * GL / 2;
						tens[a][2] = ((double)i / PP3) * GL / 2;
						a +=1;
					}
		}
		a = 0;
		if (calib_point == 189)
			for (iz = -1; iz <= 1; iz++)
				for (iy = -4; iy <= 4; iy++)
					for (ix = -3; ix <= 3; ix++) {
						tens[a][0] = (ix + (i_ch - 1) / 4 - 1) * 25.0;
						tens[a][1] = (iy + (i_ch - 1) % 4 - 1) * 25.0;
						tens[a][2] = iz * 50.0;
						a++;
					}

	////////////////////////////////////////////////////////////////////////////

	/////収束の重みの用意//////////////////////////////////////////////////////
		j = 1;
		for (i = 0; i < PATERN; i++) {
			if (i % 4 == 0)
				j++;
			alpha[i] = 1.0 / pow(10.0, j + 4.0) * (4.0 - (i % 4)) * 2.5;
			printf("%f\n", alpha[i]);
		}
	/////////////////////////////////////////////////////////////////////////////////////////
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			for (j = 0; j < 5; j++)
				dd[i_fnum][j] = 0.0;
		if (aiueo == 1)
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
				printf("Channel_No%d_File_No%d:type_theta1_theta2\n", i_ch, i_fnum + 1);
				scanf("%lf", &dd[i_fnum][3]);
				scanf("%lf", &dd[i_fnum][4]);
			}
		for (i = 0; i < NUM_SEN_COIL; i++)
			for (j = 0; j < 6; j++)
				S_p[i][j] = Scoil_pole[i][j];
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			for (i = 0; i < NUM_SEN_COIL; i++)
				r2[i_fnum][i] = 0.0;
		d_sum=0.0;
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			for (i = 0; i < NUM_SEN_COIL; i++)
				sns[i_fnum][i] = 1.0;//校正に使わないデータを0にする
			//if(i==3)
			//sns[N][i]=0.0;
		if (w_Switch == 0)
			weight(calib_buf, w);//ゲインを基に重み付け
		else
			for (i = 0; i < NUM_SEN_COIL; i++)
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					w[i_fnum][i] = 1.0;
	////////////////////////////////
	//main_loop/////////////////////////
		if (s_flg == 1)
			LOOP_MAX = 40000;//校正済み位置を使う場合はループは少なめ
		else {
			LOOP_MAX = 5000000;
			// LOOP_MAX = 1000000;
		}

		sum_buf_10000 = 1.0;
		sum_buf = 0.0;
		count1 = 0;
		count2 = 0;
		for (index = 0; index < LOOP_MAX; index++) {
			if (count1 == 5) {
			//puts("up\n");
				for (i = 0; i < PATERN; i++)
					alpha[i] *= 1.1;//収束因子を大きくする
			//printf("%f\n",alpha[i]);
				count1 = 0;
			}
			if (count2 == 5) {
			//puts("down\n");
				for (i = 0; i < PATERN; i++)
					alpha[i] /= 1.1;//収束因子を小さくする
			//printf("%f\n",alpha[i]);
				count2 = 0;
			}
			for (i = 0; i < NUM_SEN_COIL; i++) {
				for (j = 0; j < 6; j++) {
					q_a[i][j] = 0.0;
					q_b[i][j] = 0.0;
				}
				q_c[i] = 0.0;
				q_d[i] = 0.0;
			}
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
				for (i = 0; i < NUM_SEN_COIL; i++)
					for (j = 0; j < 5; j++) {
						q_a2[i_fnum][i][j] = 0.0;
						q_b2[i_fnum][i][j] = 0.0;
					}
			//各微分係数の準備////////////////////////
			for (i = 0; i < NUM_SEN_COIL; i++) {			
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {	
					for (k = 0; k < calib_point; k++) {
						if (sns[i_fnum][i] == 0.0)
							break;
						Dfuncg_dd(tens, 1.0, S_p[i], Z[i], k, i, dd[i_fnum], i_fnum * calib_point);
						//printf("Z[%d]:%lf\n",k,Z[i][k]);
						SDM_1(S_p[i], dd[i_fnum], tens, A[i][k + i_fnum * calib_point], k);
						SDM_2(tens, dd[i_fnum], k, S_p[i], A2[i_fnum][i][k + i_fnum * calib_point]);
						for (j = 0; j < 6; j++) {
							q_a[i][j] += calib_buf[i][k + i_fnum * calib_point] * A[i][k + i_fnum * calib_point][j];
							q_b[i][j] += 2 * Z[i][k + i_fnum * calib_point] * A[i][k + i_fnum * calib_point][j];
						}
						for (j = 0; j < 5; j++) {
							q_a2[i_fnum][i][j] += calib_buf[i][k + i_fnum * calib_point] * A2[i_fnum][i][k + i_fnum * calib_point][j];
							q_b2[i_fnum][i][j] += 2 * Z[i][k + i_fnum * calib_point] * A2[i_fnum][i][k + i_fnum * calib_point][j];
					//printf("q_a2[%d]:%lf\n",j,q_a2[N][i][j]);
					//printf("q_b2[%d]:%lf\n",j,q_b2[N][i][j]);
						}
						q_c[i] += pow(Z[i][k + i_fnum * calib_point], 2.0);
						q_d[i] += calib_buf[i][k + i_fnum * calib_point] * Z[i][k + i_fnum * calib_point];
					}
				}
			}
			///////////////////////////////////////////////////////////////
			for (i = 0; i < NUM_SEN_COIL; i++)
				Q[i] = q_d[i] / q_c[i];//ゲインの計算
			for (i = 0; i < NUM_SEN_COIL; i++)
				for (j = 0; j < 6; j++)
					q_dd1[i][j] = (q_a[i][j] * q_c[i] - q_d[i] * q_b[i][j]) / pow(q_c[i], 2.0);//ゲインの微分係数計算１		
			//printf("q:%lf:%lf:%lf:%lf\n",q_a[i][j],q_b[i][j],q_c[i],q_d[i]);
				//printf("q_dd1:%lf\n",q_dd1[i][j]);
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
				for (i = 0; i < NUM_SEN_COIL; i++)
					for (j = 0; j < 5; j++) {
						if (sns[i_fnum][i] == 0.0)
							break;
						q_dd2[i_fnum][i][j] = (q_a2[i_fnum][i][j] * q_c[i] - q_d[i] * q_b2[i_fnum][i][j]) / pow(q_c[i], 2.0);//ゲインの微分係数計算２
				//printf("q_dd2[%d]:%lf\n",j,q_dd2[N][i][j]);
					}
			//for(i=0;i<6;i++)printf("Q:\t%lf\n",Q[i]);
			for (i = 0; i < NUM_SEN_COIL; i++)
				for (k = 0; k < 6; k++)
					d_f[i][k] = 0.0;
			///オフセットの計算///////////////////////////
			if (rr_Switch == 0) {
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					for (i = 0; i < NUM_SEN_COIL; i++)
						r[i_fnum][i] = 0;
				if(index > 1000) {
					for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
						for (i = 0; i < NUM_SEN_COIL; i++)
							for (j = 0; j < calib_point; j++)
								r[i_fnum][i] += calib_buf[i][j + i_fnum * calib_point] - Q[i]*Z[i][j + i_fnum * calib_point];
					for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
						for (i = 0; i < NUM_SEN_COIL; i++)
							r[i_fnum][i] = r[i_fnum][i] / (double)((calib_point));			
					if(r_Switch == 0)
						for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
							for (i = 0; i < NUM_SEN_COIL; i++) {
								for (j = 0; j < calib_point; j++)
									calib_buf[i][j + i_fnum * calib_point] -= r[i_fnum][i];
								r2[i_fnum][i] += r[i_fnum][i];
							}
				}
			} else {
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					for (i = 0; i < NUM_SEN_COIL; i++)
						r[i_fnum][i] = 0;
				if (index > 1000) {
					for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
						for (i = 0; i < NUM_SEN_COIL; i++)
							for (j = 0; j < calib_point; j++)
								r[0][i] += calib_buf[i][j + i_fnum * calib_point] - Q[i]*Z[i][j + i_fnum * calib_point];
					for (i = 0; i < NUM_SEN_COIL; i++)
						r[0][i] = r[0][i] / (double)(calib_point * FILENUM);
					//	printf("p:%f\n",r[i]);
					if (r_Switch == 0)
						for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
							for (i = 0; i < NUM_SEN_COIL; i++) {
								for (j = 0; j < calib_point; j++)
									calib_buf[i][j + i_fnum * calib_point] -= r[0][i];
								r2[i_fnum][i] += r[0][i];
							}
				}
			}			
		/////////////////////////////////////////////////////////////////////////////////////		
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
				for (j = 0; j < 5; j++)
					d_f3[i_fnum][j] = 0.0;
				//printf("d_f3[%d]:%lf:%lf\n",N,d_f3[N][3],d_f3[N][4]);
				for (i = 0; i < NUM_SEN_COIL; i++) {
					for (k = 0; k < calib_point; k++) {		
						if (sns[i_fnum][i] == 0.0)
							break;
						for (j = 0; j < 6; j++)
							if (index > 1999)
								d_f[i][j] += w[i_fnum][i] * sns[i_fnum][i] * -2.0 * (calib_buf[i][k + i_fnum * calib_point] - Q[i] * Z[i][k + i_fnum * calib_point]) * (q_dd1[i][j] * Z[i][k + i_fnum * calib_point] + Q[i] * A[i][k + i_fnum * calib_point][j]);
						//送信コイル位置位置の微分係数
						if (i_fnum == 1000)
							for (j = 3; j < 5; j++)
								d_f3[i_fnum][j] += w[i_fnum][i] * sns[i_fnum][i] * -2.0 * (calib_buf[i][k + i_fnum * calib_point] - Q[i] * Z[i][k + i_fnum * calib_point]) * (q_dd2[i_fnum][i][j] * Z[i][k + i_fnum * calib_point] + Q[i] * A2[i_fnum][i][k + i_fnum * calib_point][j]);
						//受信コイル角度の微分係数	
						else
							for (j = 0; j < 5; j++)
								d_f3[i_fnum][j] += w[i_fnum][i] * sns[i_fnum][i] * -2.0 * (calib_buf[i][k + i_fnum * calib_point] - Q[i] * Z[i][k + i_fnum * calib_point]) * (q_dd2[i_fnum][i][j] * Z[i][k + i_fnum * calib_point] + Q[i] * A2[i_fnum][i][k + i_fnum * calib_point][j]);
						//受信コイル位置と角度の微分係数
					}
				}
			}
			if (Switch == 1)
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					for (j = 0; j < 3; j++)
						d_f3[i_fnum][j] = 0.0;
			if (s_flg == 1)
				for (i = 0; i < NUM_SEN_COIL; i++)
					for (j = 0; j < 6; j++)
						d_f[i][j] = 0.0;
			////////////////////////////////////////////////////////////////////
			//////収束因子ごとの更新後誤差計算////////////////////////////////
			for (k = 0; k < PATERN; k++) {
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
					if (i_fnum == 1000) {
						for (j = 0; j < 3; j++)
							dd_buf[k][i_fnum][j] = 0.0;
						for (j = 3; j < 5; j++)
							dd_buf[k][i_fnum][j] = dd[i_fnum][j] - alpha[k] * d_f3[i_fnum][j];
					} else {
						for (j = 0; j < 5; j++)
							dd_buf[k][i_fnum][j] = dd[i_fnum][j] - alpha[k] * d_f3[i_fnum][j];
					}
				}
				for (i = 0; i < NUM_SEN_COIL; i++)
					for (j = 0; j < 6; j++)
						S_p_buf[i][k][j] = S_p[i][j] - alpha[k] * d_f[i][j];
				sum2[k] = 0.0;
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					for (i = 0; i < NUM_SEN_COIL; i++)
						for (j = 0; j < calib_point; j++)
							Dfuncg_dd(tens, 1.0, S_p_buf[i][k], Z[i], j, i, dd_buf[k][i_fnum], i_fnum * calib_point);
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					for (i = 0; i < NUM_SEN_COIL; i++) {
						sum3[i][i_fnum][k] = 0.0;
						for (j = 0; j < calib_point; j++) {
							if (sns[i_fnum][i] == 0.0)
								break;
							sum2[k] += pow((calib_buf[i][j + (i_fnum * calib_point)] - Q[i] * Z[i][j + (i_fnum * calib_point)]), 2.0);//全データの残差誤差
							sum3[i][i_fnum][k] += pow((calib_buf[i][j + (i_fnum * calib_point)]- Q[i] * Z[i][j + (i_fnum * calib_point)]), 2.0);//送信コイル・データセットごとの残差誤差
						}
					}
			}
			MIN = min_index(sum2, 0, PATERN);//収束因子の決定
			if (MIN == 0)
				count1++;
			else
				count1 = 0;
			if (MIN == 7)
				count2++;
			else
				count2 = 0;
			//受信コイルの位置と角度の更新/////////////////////////////////////////
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
				if (i_fnum == 1000)
					for (j = 3; j < 5; j++)
						dd[i_fnum][j] = dd[i_fnum][j] - alpha[MIN] * d_f3[i_fnum][j];
				else
					for (j = 0; j < 5; j++)
						dd[i_fnum][j] = dd[i_fnum][j] - alpha[MIN] * d_f3[i_fnum][j];			
			}
		//////////////////////////////////////////////////////////////////////////////////
			for (i = 0; i < NUM_SEN_COIL; i++)
				for (j = 0; j < 6; j++)
					S_p[i][j] = S_p[i][j] - alpha[MIN] * d_f[i][j];//送信コイル位置の更新
			d_sum = 0.0;
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
				for (j = 0; j < 5; j++)
					d_sum += fabs(d_f3[i_fnum][j]);	
			for (i = 0; i < NUM_SEN_COIL; i++)
				for (j = 0; j < 6; j++)
					d_sum += fabs(d_f[i][j]);
			///各変数の表示
			if (index % 100 == 0) {
				for (i = 0; i < NUM_SEN_COIL; i++)
					printf("d_f[%d]:%lf:%lf:%lf:%lf:%lf:%lf\n", i, d_f[i][0], d_f[i][1], d_f[i][2], d_f[i][3], d_f[i][4], d_f[i][5]);
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					printf("d_f3[%d]:%lf:%lf:%lf:%lf:%lf\n", i_fnum, d_f3[i_fnum][0], d_f3[i_fnum][1], d_f3[i_fnum][2], d_f3[i_fnum][3], d_f3[i_fnum][4]);
				for (i = 0; i < NUM_SEN_COIL; i++) {
					printf("S_pp[%d]:%lf:%lf;%lf\n", i, S_p[i][0], S_p[i][1], S_p[i][2]);
					printf("S_pn[%d]:%lf;%lf;%lf\n", i, S_p[i][3], S_p[i][4], S_p[i][5]);
				}
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
					printf("d_f:%lf:%lf:%lf:%lf:%lf\n", dd[i_fnum][0], dd[i_fnum][1], dd[i_fnum][2], dd[i_fnum][3], dd[i_fnum][4]);
				printf("sum2[%d]\t%f\n", index, sum2[MIN]);
				printf("d_sum:%lf\n", d_sum / (FILENUM * 5.0 + 36.0));
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
					printf("[%d][%d]", i_fnum, i);
					for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
						printf(":%lf", sum3[i_send][i_fnum][MIN]);
					printf("\n");
					sum_sum3 = 0.0;
					for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
						sum_sum3 += sum3[i_send][i_fnum][MIN];
					printf("[%d]:%lf\n", i_fnum, sum_sum3);
				}
				printf("MIN:%d:%lf\n", MIN, alpha[MIN]);
				for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
					printf("r2");
					for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
						printf(":%lf", r2[i_fnum][i_send]);
					printf("\n");
				}
				printf("test:%lf\n", sum2[MIN] / sum_buf_10000);
			}
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
				sum4[i_fnum] = 0.0;
				for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
					sum4[i_fnum] += sum3[i_send][i_fnum][MIN];
			}
			if (index == 0)
				sum_buf_10000 = sum2[MIN];
			if (index > 999)
				if (index % 1000 == 0) {
					if((sum2[MIN] / sum_buf_10000) > 0.99999) {
						puts("LOOOPEND\n");
						break;
					}
					sum_buf_10000 = sum2[MIN];
				}
			if (index > 10000 && alpha[MIN] < 1.0 / 1000000000.0) {
				puts("LOOOPEND2\n");
				break;
			}
			sum_buf = sum2[MIN];
			for (i_fnum = 0; i_fnum < FILENUM; i_fnum++) {
				for (i = 0; i < NUM_SEN_COIL; i++)
					sum3_buf[i][i_fnum] = sum3[i][i_fnum][MIN];
				sum4_buf[i_fnum] = 0.0;
				for (i_send = 0; i_send < NUM_SEN_COIL; i_send++)
					sum4_buf[i_fnum] += sum3[i_send][i_fnum][MIN];//送信コイルごとのデータセット
			}
		}
		//////////////////////////////////
		//file_write//
		fprintf(fp_res, "sum_buf:\t%f\n", sum_buf);
		for (i = 0; i < NUM_SEN_COIL; i++)
			fprintf(fp_res, "Q:\t%f\n", Q[i]);
		for (i = 0; i < NUM_SEN_COIL; i++) {
			for (j = 0; j < 3; j++)
				S_c[i][j] = (S_p[i][j] + S_p[i][j + 3]) / 2.0;
			//printf("S_c[%d][%d]:%f\n",i,j,S_c[i][j]);
			S_l[i] = pow((pow(S_p[i][0] - S_p[i][3], 2.0) + pow(S_p[i][1] - S_p[i][4], 2.0) + pow(S_p[i][2] - S_p[i][5], 2.0)), 0.5);
			printf("S_l[%d]%lf\n", i, S_l[i]);
			fprintf(fp_res, "S_l[%d]\t%lf\n", i, S_l[i]);
			for (j = 0; j < 3; j++) {
				S_b[i][j] = S_p[i][j] - S_c[i][j];
				S_b[i][j] = S_b[i][j] / (S_l[i] / 2.0);
				printf("S_b[%d][%d]:%lf\n", i, j, S_b[i][j]);
			}
			theta[i][0] = -asin(S_b[i][2]);
			theta[i][1] = acos(S_b[i][0] / cos(theta[i][0]));
			printf("theta[%d]\t%lf\t%lf\n", i, theta[i][0] * 180.0 / M_PI, theta[i][1] * 180.0 / M_PI);		 
			fprintf(fp_res, "theta[%d]\t%lf\t%lf\n", i, theta[i][0] * 180.0 / M_PI, theta[i][1] * 180.0 / M_PI);	
		}
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			fprintf(fp_res, "Scoil_res\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", sum3[0][i_fnum][MIN], sum3[1][i_fnum][MIN], sum3[2][i_fnum][MIN], sum3[3][i_fnum][MIN], sum3[4][i_fnum][MIN], sum3[5][i_fnum][MIN], sum3[6][i_fnum][MIN], sum3[7][i_fnum][MIN]);
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			fprintf(fp_res, "r2\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", r2[i_fnum][0], r2[i_fnum][1], r2[i_fnum][2], r2[i_fnum][3], r2[i_fnum][4], r2[i_fnum][5], r2[i_fnum][6], r2[i_fnum][7]);
		for (i = 0; i < NUM_SEN_COIL; i++)
			fprintf(fp_Sp, "%f\n", Q[i]);
		for (i = 0; i < NUM_SEN_COIL; i++)
			fprintf(fp_Sp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", S_p[i][0], S_p[i][1], S_p[i][2], S_p[i][3], S_p[i][4], S_p[i][5]);
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			fprintf(fp_Sp, "%lf\t%lf\t%lf\t%lf\t%lf\n", dd[i_fnum][0], dd[i_fnum][1], dd[i_fnum][2], dd[i_fnum][3], dd[i_fnum][4]);
		fprintf(fp_Sp, "%d\n", index);
		for (i = 0; i < NUM_SEN_COIL; i++)
			for (j = 0; j < 2; j++)
				if (theta[i][j] > 0)
					theta2[i][j] = theta[i][j] + M_PI / 2.0;
				else
					theta2[i][j] = theta[i][j] + M_PI / 2.0;
		for (i = 0; i < NUM_SEN_COIL; i++)
			fprintf(fp_Sc, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", S_p[i][0], S_p[i][1], S_p[i][2], S_p[i][3], S_p[i][4], S_p[i][5]);
		fprintf(fp_GP, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", Q[0], Q[1], Q[2], Q[3], Q[4], Q[5], Q[6], Q[7]);
		fprintf(fp_GP, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", r2[0][0], r2[0][1], r2[0][2], r2[0][3], r2[0][4], r2[0][5], r2[0][6], r2[0][7]);
		for (i_fnum = 0; i_fnum < FILENUM; i_fnum++)
			if(in_fp[i_fnum] != NULL)
				fclose(in_fp[i_fnum]);
		//exit(0);
		if (fp_res != NULL)
			fclose(fp_res);
		if (fp_Sp != NULL)
			fclose(fp_Sp);
		if (fp_GP != NULL)
			fclose(fp_GP);
		if (fp_Sc != NULL)
			fclose(fp_Sc);
		if(s_flg ==1)
			if(fp_Sc_in != NULL)
				fclose(fp_Sc_in);
	}
	/////////////////////////////////
	printf("program is ending...type some number");
	scanf("%d",&i);

	free_matrix(calib_buf,NUM_SEN_COIL);
	free_matrix(Z,NUM_SEN_COIL);
	free_matrix(Z2,NUM_SEN_COIL);
	free_matrix(tens,calib_point);
	free_matrix(Scoil_pole,NUM_SEN_COIL);
	free_matrix(S_p,NUM_SEN_COIL);
	free_matrix(w,FILENUM);


	//exit(0);


	return (0);

}