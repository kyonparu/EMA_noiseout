// sumi2.cpp : �R���\�[�� �A�v���P�[�V�����̃G���g�� �|�C���g���`���܂��B
//

/*
#define _USE_MATH_DEFINES
*/

#include "stdafx.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<float.h>

#define CL 60.0    //�R�C����[m]
#define L  337.5    //���W���S����R�C���̒��S�܂ł̒���[m]
#define K 8         // ���M�`���l����
/*
#define M_PI 3.141592653589793
*/
#define PATERN 1
#define FILENUM 1  //�f�[�^�Z�b�g�̐�
//#define calib_point 27
#define calib_point 189
#define Switch 0
// #define N 27//�f�[�^�Z�b�g������̃f�[�^��
#define v_switch 0//�X���̐؂�ւ��@����0�\��
#define START_CH_NO 1 //�m�C�Y�����`���l��No.�̍ŏ��l
#define END_CH_NO 12 //�m�C�Y�����`���l��No.�̍ő�l
#define REF_CH_NO 2 //��ʑ��v�Z�p�`���l��No.



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
void sai1(double *x,double *y,double *z,double *z2,int I)
{
        int i;
        double a0,a1,p,q;
        double A00,A01,A02,A11,A12;

        //FILE *output1;
       // FILE *output2;

        //output1=fopen("output1.data","w");
        //output2=fopen("output2.data","w");


		
	//puts("ok");

	//scanf("%d",&i);



        A00=A01=A02=A11=A12=0.0;

        for (i=0;i<calib_point;i++) {
                A00+=1.0;
                A01+=x[i];
                A02+=y[i];
                A11+=x[i]*x[i];
                A12+=x[i]*y[i];
        }

		/*�P�����̌W���̌v�Z*/
        a0=(A02*A11-A01*A12)/(A00*A11-A01*A01);
        a1=(A00*A12-A01*A02)/(A00*A11-A01*A01);

		z[I]=a1;
		z2[I]=a0;

/*gnuplot�ŃO���t�\���̂��߂Ƀf�[�^��t�@�C���ɕۑ�*/
       // for(i=0;i<N;i++) {
        //        fprintf(output1,"%f %f\n",x[i],y[i]);
     //   }

/*gnuplot�ŃO���t�\���̂��߂ɁA�v�Z�œ���ꂽ�P�����Ńf�[�^�ߕӂ̃O���t�f�[�^��ۑ�*/
        for(p=x[0]-10.0;p<x[calib_point-1]+10.0;p+=0.01) {
                q=a0+a1*p;
         //       fprintf(output2,"%f %f\n",p,q);
        }

       // fclose(output1);
       // fclose(output2);

}
void sort(double **new_d_buf,double **d_buf)
{
	int i,j,k;
	int MIN;

	for(i=0; i < K; i++)
	{
		for(j=0;j<calib_point;j++)
		{

			MIN=min_index(d_buf[2*i],0,calib_point);

			new_d_buf[2*i][j]=d_buf[2*i][MIN];
			new_d_buf[2*i+1][j]=d_buf[2*i+1][MIN];

			d_buf[2*i][MIN]=DBL_MAX;
			//printf("%f:%f\n",new_d_buf[2*i][j],new_d_buf[2*i+1][j]);
		}
	}

}
void fitting(double **cmp_data,double *s1,double *s2)
{
	int i,j,k;

	for(i=0; i < K; i++)
	{
		for(j=0;j<calib_point;j++)
		{
			cmp_data[2*i][j]=cmp_data[2*i][j] - s1[i];
			cmp_data[2*i+1][j]=cmp_data[2*i+1][j] - s2[i];
		}
	}
}

void fitting2(double **cmp_data,double *v,double **mean_data)
{
	int i,j,k,s;

	for(s=0;s<FILENUM;s++)
	{
		for(i=0; i < K; i++)
		{
			for(j=0;j<calib_point;j++)
			{		
				mean_data[i][j+s*calib_point] = (cmp_data[2*i][j+s*calib_point] + cmp_data[2*i+1][j+s*calib_point]*v[i])
	           		                            /pow((1.0+pow(v[i],2.0)),0.5);
			}
		}
	}
}

void fitting3(double **cmp_data,double **cst,double **mean_data)
{
	int i,j,k,s;

	for(s=0;s<FILENUM;s++)
	{
		for(i=0; i < K;i++)
		{
			for(j=0;j<calib_point;j++)
			{
				mean_data[i][j+s*calib_point]=(cmp_data[2*i][j+s*calib_point] + cmp_data[2*i+1][j+s*calib_point]*cst[s][i])/pow((1.0+pow(cst[s][i],2.0)),0.5);
			}
		}
	}
}


int main()
{


	FILE *in_fp[FILENUM],*out_fp[FILENUM],*old_fp[FILENUM],*cst_fp;

	int i, j, s, i_ch, i_send;
	int use_cst_file, filenum_in_cst, use_file_no;
	char di_dir[FILENUM][256];
	char d_name[256], di_name[256], cst_dir_name[256], cst_file_name[256];
	char f_name[256];
	double **d_buf,**d_buf2;
	double **new_d_buf;
	double **cmp_data,**mean_data,**old_mean_data;
	double **cst,**cst2;
	double v[K];
	double sum=0;

	/* �̈�m�� */
	mean_data     =  new_matrix(K,calib_point*FILENUM);
	old_mean_data =  new_matrix(K,calib_point*FILENUM);
	d_buf         =  new_matrix(2 * K,calib_point);
	d_buf2        =  new_matrix(2 * K,calib_point);
	new_d_buf     =  new_matrix(2 * K,calib_point);
	cmp_data      =  new_matrix(2 * K,calib_point*FILENUM);
	cst           =  new_matrix(FILENUM, K);
	cst2          =  new_matrix(FILENUM, K);

	for (i = 0; i < FILENUM; i++) {
		printf("put_input_file_dir(No. %d)\n", i+1);
		scanf("%s", di_dir[i]);
	}

	printf("Do you use the existed cst file?(Yes = 1, No = others)\n");
	scanf("%d", &use_cst_file);
	if (use_cst_file == 1) {
		printf("Put the number of FILES for calculation in cst file?\n");
		scanf("%d", &filenum_in_cst);
		printf("Which number of FILE is used for noise reduction\n");
		scanf("%d", &use_file_no);
		printf("put_cst_dir_name\n");
		scanf("%s", cst_dir_name);
	}

	for (i_ch = START_CH_NO; i_ch <= END_CH_NO; i_ch++) {
		if (i_ch == REF_CH_NO)
			continue;
		for (s = 0; s < FILENUM; s++) {
			sprintf(di_name, "..\\..\\data\\%s\\c-mean\\ch%d.data", di_dir[s], i_ch);
			sprintf(d_name, "..\\..\\data\\%s\\c-mean\\out_ch%d", di_dir[s], i_ch);
			in_fp[s] = fopen(di_name, "r");
			sprintf(f_name, "%s.data", d_name);
			out_fp[s] = fopen(f_name, "w");
			sprintf(f_name, "%s-old.data", d_name);
			old_fp[s] = fopen(f_name, "w");
		}
		if (use_cst_file == 1) {
			sprintf(cst_file_name, "..\\..\\data\\%s\\c-mean\\out_ch%d", cst_dir_name, i_ch);
			sprintf(f_name, "%s_cst.data", cst_file_name);
			cst_fp = fopen(f_name, "r");
			for (s = 0; s < filenum_in_cst; s++) {
				for (i = 0; i < K; i++)
					fscanf(cst_fp, "%lf\t", &cst[s][i]);
				fscanf(cst_fp, "\n");
			}
			fscanf(cst_fp,"\n");
			for (s = 0; s< filenum_in_cst; s++) {
				for (i = 0; i < K; i++)
					fscanf(cst_fp, "%lf\t", &cst2[s][i]);
				fscanf(cst_fp,"\n");
			}
		} else {
			sprintf(f_name,"%s_cst.data",d_name);
			cst_fp=fopen(f_name,"w");
		}
		for (s = 0; s < FILENUM; s++) {
		/*(����+����)*K �̐��l��buf�Ɋi�[ */
			for (i = 0; i < calib_point; i++) {
				for (i_send = 0; i_send < K - 1; i_send++)
					fscanf(in_fp[s], "%lf\t%lf\t", &d_buf[2 * i_send][i], &d_buf[2 * i_send + 1][i]);
				fscanf(in_fp[s], "%lf\t%lf\n", &d_buf[2 * (K - 1)][i], &d_buf[2 * (K - 1) + 1][i]);
				for (j = 0; j < 2 * K; j++)
					d_buf2[j][i] = d_buf[j][i];
			}
			if (use_cst_file != 1) {
				/* ��A������������߂Ɂi�����A�����j�̃f�[�^������̑傫�����ɕ��ёւ� */
				sort(new_d_buf, d_buf);
				for (i = 0; i < K; i++) {
					sai1(new_d_buf[2 * i], new_d_buf[2 * i + 1], cst[s], cst2[s], i); //��A�����̌X���ƐؕЂ�v�Z
					printf("%lf\t%lf\n", cst[s][i], cst2[s][i]);
				//s1[i] = (-1.0*cst[i][0]*cst[i][1]/(cst[i][0]*cst[i][0]+1.0));
				///s2[i] = (cst[i][1]/(cst[i][0]*cst[i][0]+1.0));
				}
			}
			for (i = 0; i < 2 * K; i++)
				for (j = 0; j < calib_point; j++) 
					cmp_data[i][j + s * calib_point] = d_buf2[i][j];
				//printf("%lf\n",d_buf2[i][j]);
		///exit(0);
		}
		//exit(0);
		for (s = 0; s < FILENUM; s++)
			for (j = 0; j < K; j++)
				for (i = 0; i < calib_point; i++) {
					old_mean_data[j][i + s * calib_point] = pow((pow(cmp_data[2 * j][i + s * calib_point], 2.0) 
						+ pow(cmp_data[2 * j + 1][i + s * calib_point], 2.0)), 0.5);//�]���̎�M�M���l�̌v�Z
					if(cmp_data[2 * j][i + s * calib_point] < 0.0)
						old_mean_data[j][i + s * calib_point] *= -1.0;
				}
				//exit(0);
		if (use_cst_file != 1) {
			for (i = 0; i < K; i++) {
				sum = 0.0;
				for (s = 0; s < FILENUM; s++)
					sum += cst[s][i];
				v[i]=sum/(double)FILENUM;
				//printf("V:%f\n",v[i]);
			}
			if (v_switch == 0)
				fitting2(cmp_data, v, mean_data);//�f�[�^�Z�b�g�Ԃ̕��ς̌X���Ŏʑ�
			else
				fitting3(cmp_data, cst, mean_data);//�f�[�^�Z�b�g���̌X���Ŏʑ�
		} else
			fitting2(cmp_data, cst[use_file_no - 1], mean_data);
		//exit(0);
		//�t�@�C���o��///
		for (s = 0; s < FILENUM; s++) {
			for (i = 0; i < calib_point; i++) {
				for (i_send = 0; i_send < K - 1; i_send++)
					fprintf(out_fp[s], "%lf\t", mean_data[i_send][i + s * calib_point]);
				fprintf(out_fp[s], "%lf\n", mean_data[K - 1][i + s * calib_point]);
				for (i_send = 0; i_send < K - 1; i_send++)
					fprintf(old_fp[s], "%lf\t", old_mean_data[i_send][i + s * calib_point]);
				fprintf(old_fp[s], "%lf\n", old_mean_data[K - 1][i + s * calib_point]);
			}
			fprintf(out_fp[s], "%lf\t%lf\n", 0.0, 0.0);
			fprintf(old_fp[s], "%lf\t%lf\n", 0.0, 0.0);
		}
		if (use_cst_file != 1) {
			for (s = 0; s < FILENUM; s++) {
				for (i = 0; i < K; i++)
					fprintf(cst_fp, "%lf\t", cst[s][i]);
				fprintf(cst_fp,"\n");
			}
			fprintf(cst_fp,"\n");
			for (s = 0; s < FILENUM; s++) {
				for (i = 0; i < K; i++)
					fprintf(cst_fp, "%lf\t", cst2[s][i]);
				fprintf(cst_fp,"\n");
			}
		}
		//exit(0);
		//scanf("%d",&i);
		for (s = 0; s < FILENUM; s++) {
			fclose(in_fp[s]);
			fclose(out_fp[s]);
			fclose(old_fp[s]);
		}
		fclose(cst_fp);
	}

	free_matrix(d_buf, 2 * K);	
	free_matrix(d_buf2, 2 * K);
	free_matrix(new_d_buf, 2 * K);
	free_matrix(mean_data, K);
	free_matrix(old_mean_data, K);
	free_matrix(cst,FILENUM);
	free_matrix(cmp_data, 2 * K);
	free_matrix(cst2,FILENUM);

        return 0;
}