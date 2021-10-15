#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <x86intrin.h>
#include <stdint.h>


/*************************************************
* Name:         B1-Hybrid1 Algorithm corresponds to b1_real() below.
*               See the paper: "Improved polynomial multiplication algorithms over characteristic three fields and applications to NTRU Prime" 
                                by E. Yeniaras and M. Cenk for the details of the algorithm:
*              
*
* Description: Multiplies two polynomials with input size "n=654" with coefficients in F_3
*              Using 1 level of Bernstein's 3-way algorithm B1, then Karatsuba 2-way and the unbalanced Karatsuba 2-way down to the input sizes 16 
               and then the schoolbook at the final level.
*
* Arguments:   struct complex a[]:        array of F_3 coefficients for the input polynomial a(x)
*              struct complex b[]:        array of F_3  coefficients for the input polynomial b(x)
*              struct complex c[]:        array of F_3  coefficients for the multiplication polynomial c(x)=a(x).b(x)    
*              int n=654                  number of coefficients for input polynomials a(x) and b(x) 
**************************************************/



struct complex{

	int re;
	int im;
};

void copy_same_loc(struct complex *source, struct complex *destination, int from, int to){

	for(int i=from; i<to; ++i)
		destination[i] = source[i];
}

void copy(struct complex *source, struct complex *destination, int from, int to){

	for(int i=from, j=0; i<to; ++i, ++j)
		destination[j] = source[i];
}

void sum(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].re;
		target[i].im = a[i].im + b[i].im;
	}
}

void sum_neg_neg(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].re - b[i].re;
		target[i].im = - a[i].im - b[i].im;
	}
}

void sum_multiply_w(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].im - b[i].im;
		target[i].im = a[i].re + b[i].re;
	}
}

void difference(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].re;
		target[i].im = a[i].im - b[i].im;
	}
}

void multiplication(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		int ac = a[i].re * b[i].re;
		int bd = a[i].im * b[i].im;
		target[i].re = ac - bd;
		target[i].im = (a[i].re + a[i].im)*(b[i].re + b[i].im) - ac - bd;
	}
}

void negation(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = -a[i].re;
		target[i].im = -a[i].im;
	}
}


void sum_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].re;

	}
}

void sum_neg_neg_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = - a[i].re - b[i].re;

	}
}

void sum_re_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].im;

	}
}

void sum_re_neg_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].im;

	}
}

void sum_im_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].im + b[i].im;

	}
}

void sum_im_neg_im_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].im - b[i].im;

	}
}

void difference_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].re;

	}
}

void multiplication_real(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re * b[i].re;

	}
}

void negation_real(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = -a[i].re;

	}
}

void create_term(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re - b[i].im;
		target[i].im = a[i].im + b[i].re;
	}
}

void create_term_neg(struct complex *a, struct complex *b, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re + b[i].im;
		target[i].im = a[i].im - b[i].re;
	}
}

void phase_change(struct complex *a, struct complex *target, int n){

	for(int i=0; i<n; ++i){
		target[i].re = a[i].re;
		target[i].im = -a[i].im;		
	}
}

void sb_comba_complex(struct complex a[], struct complex b[], struct complex c[], int n){
	
	int re;
	int im;

	for(int i=0; i<n; ++i){
		re = im = 0;
		for(int j=0; j<=i; ++j){
			re += a[j].re * b[i-j].re - a[j].im * b[i-j].im;
			im += a[j].re * b[i-j].im + a[j].im * b[i-j].re;
		}
		c[i].re = re;
		c[i].im = im;
	}

	for(int i=n; i<2*n-1; ++i){
		re = im = 0;
		for(int j=i-n+1; j<n; ++j){
			re += a[j].re * b[i-j].re - a[j].im * b[i-j].im;
			im += a[j].re * b[i-j].im + a[j].im * b[i-j].re;
		}
		c[i].re = re;
		c[i].im = im;
	}

}

void sb_comba_real(struct complex a[], struct complex b[], struct complex c[], int n){

	int result;
int i,j;
	for( i=0; i<n; ++i){
		result = 0;
		for( j=0; j<=i; ++j)
			result += a[j].re*b[i-j].re;
		c[i].re = result;
	}

	for( i=n; i<2*n-1; ++i){
		result = 0;
		for( j=i-n+1; j<n; ++j)
			result += a[j].re*b[i-j].re;
		c[i].re = result;
	}

}



void ub_real_27(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define UB_REAL_MAX_AB_SIZE 14 // ceil(n/2)
	#define UB_REAL_MAX_P_SIZE 27 // n

	/* base case */
	if(n == 1){
		multiplication_real(a, b, c, n);
		return;
	}

	int n0 = (n+1)/2;
	int n1 = n-n0;

	struct complex *a0 = a;
	struct complex *a1 = a+n0;
	struct complex *b0 = b;
	struct complex *b1 = b+n0;

	struct complex sum_a0_a1[UB_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[UB_REAL_MAX_AB_SIZE];

	// assumed n is odd.
	sum(a0, a1, sum_a0_a1, n1);
	sum_a0_a1[n0-1] = a0[n0-1];
	sum(b0, b1, sum_b0_b1, n1);
	sum_b0_b1[n0-1] = b0[n0-1];


	struct complex p1[UB_REAL_MAX_P_SIZE];
	struct complex p2[UB_REAL_MAX_P_SIZE];
	struct complex p3[UB_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	sb_comba_real(a0, b0, p1, n0);
	sb_comba_real(sum_a0_a1, sum_b0_b1, p2, n0);
	sb_comba_real(a1, b1, p3, n1);


	sum_real(c, p1, c, 2*n0-1);
	difference_real(c+n0, p3, c+n0, 2*n1-1);

	for(int i=n+n1-2; i>=0 ; --i){
		c[i+n0].re -= c[i].re;
		c[i+n0].im -= c[i].im;
	}

	sum_real(c+n0, p2, c+n0, 2*n0-1);
	
}



void ka2_real_28(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 14 // n/2
	#define KA2_REAL_MAX_P_SIZE 27 // n-1



	int n2 = n/2;

	struct complex *a0 = a;
	struct complex *a1 = a+n2;
	struct complex *b0 = b;
	struct complex *b1 = b+n2;

	struct complex sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

	sum(a0, a1, sum_a0_a1, n2);
	sum(b0, b1, sum_b0_b1, n2);

	struct complex p1[KA2_REAL_MAX_P_SIZE];
	struct complex p2[KA2_REAL_MAX_P_SIZE];
	struct complex p3[KA2_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	sb_comba_real(a0, b0, p1, n2);
	sb_comba_real(sum_a0_a1, sum_b0_b1, p2, n2);
	sb_comba_real(a1, b1, p3, n2);


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}

void ub_real_55(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define UB_REAL_MAX_AB_SIZE 28 // ceil(n/2)
	#define UB_REAL_MAX_P_SIZE 55 // n

	/* base case */
	if(n == 1){
		multiplication_real(a, b, c, n);
		return;
	}

	int n0 = (n+1)/2;
	int n1 = n-n0;

	struct complex *a0 = a;
	struct complex *a1 = a+n0;
	struct complex *b0 = b;
	struct complex *b1 = b+n0;

	struct complex sum_a0_a1[UB_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[UB_REAL_MAX_AB_SIZE];

	// assumed n is odd.
	sum(a0, a1, sum_a0_a1, n1);
	sum_a0_a1[n0-1] = a0[n0-1];
	sum(b0, b1, sum_b0_b1, n1);
	sum_b0_b1[n0-1] = b0[n0-1];


	struct complex p1[UB_REAL_MAX_P_SIZE];
	struct complex p2[UB_REAL_MAX_P_SIZE];
	struct complex p3[UB_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ka2_real_28(a0, b0, p1, n0);
	ka2_real_28(sum_a0_a1, sum_b0_b1, p2, n0);
	ub_real_27(a1, b1, p3, n1);


	sum_real(c, p1, c, 2*n0-1);
	difference_real(c+n0, p3, c+n0, 2*n1-1);

	for(int i=n+n1-2; i>=0 ; --i){
		c[i+n0].re -= c[i].re;
		c[i+n0].im -= c[i].im;
	}

	sum_real(c+n0, p2, c+n0, 2*n0-1);
	
}



void ka2_real_54(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 27 // n/2
	#define KA2_REAL_MAX_P_SIZE 53// n-1


	int n2 = n/2;

	struct complex *a0 = a;
	struct complex *a1 = a+n2;
	struct complex *b0 = b;
	struct complex *b1 = b+n2;

	struct complex sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

	sum(a0, a1, sum_a0_a1, n2);
	sum(b0, b1, sum_b0_b1, n2);

	struct complex p1[KA2_REAL_MAX_P_SIZE];
	struct complex p2[KA2_REAL_MAX_P_SIZE];
	struct complex p3[KA2_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ub_real_27(a0, b0, p1, n2);
	ub_real_27(sum_a0_a1, sum_b0_b1, p2, n2);
	ub_real_27(a1, b1, p3, n2);


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}

void ub_real_109(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define UB_REAL_MAX_AB_SIZE 55 // ceil(n/2)
	#define UB_REAL_MAX_P_SIZE 109 // n

	/* base case */
	if(n == 1){
		multiplication_real(a, b, c, n);
		return;
	}

	int n0 = (n+1)/2;
	int n1 = n-n0;

	struct complex *a0 = a;
	struct complex *a1 = a+n0;
	struct complex *b0 = b;
	struct complex *b1 = b+n0;

	struct complex sum_a0_a1[UB_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[UB_REAL_MAX_AB_SIZE];

	// assumed n is odd.
	sum(a0, a1, sum_a0_a1, n1);
	sum_a0_a1[n0-1] = a0[n0-1];
	sum(b0, b1, sum_b0_b1, n1);
	sum_b0_b1[n0-1] = b0[n0-1];


	struct complex p1[UB_REAL_MAX_P_SIZE];
	struct complex p2[UB_REAL_MAX_P_SIZE];
	struct complex p3[UB_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ub_real_55(a0, b0, p1, n0);
	ub_real_55(sum_a0_a1, sum_b0_b1, p2, n0);
	ka2_real_54(a1, b1, p3, n1);


	sum_real(c, p1, c, 2*n0-1);
	difference_real(c+n0, p3, c+n0, 2*n1-1);

	for(int i=n+n1-2; i>=0 ; --i){
		c[i+n0].re -= c[i].re;
		c[i+n0].im -= c[i].im;
	}

	sum_real(c+n0, p2, c+n0, 2*n0-1);
	
}





void ka2_real_218(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 109 // n/2
	#define KA2_REAL_MAX_P_SIZE 217 // n-1



	int n2 = n/2;

	struct complex *a0 = a;
	struct complex *a1 = a+n2;
	struct complex *b0 = b;
	struct complex *b1 = b+n2;

	struct complex sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

	sum(a0, a1, sum_a0_a1, n2);
	sum(b0, b1, sum_b0_b1, n2);

	struct complex p1[KA2_REAL_MAX_P_SIZE];
	struct complex p2[KA2_REAL_MAX_P_SIZE];
	struct complex p3[KA2_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ub_real_109(a0, b0, p1, n2);
	ub_real_109(sum_a0_a1, sum_b0_b1, p2, n2);
	ub_real_109(a1, b1, p3, n2);


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}





void ka2_real_220(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define KA2_REAL_MAX_AB_SIZE 110 // n/2
	#define KA2_REAL_MAX_P_SIZE 219 // n-1

	/* base case */
	if(n == 55){
		ub_real_55(a, b, c, n);
		return;
	}

	int n2 = n/2;

	struct complex *a0 = a;
	struct complex *a1 = a+n2;
	struct complex *b0 = b;
	struct complex *b1 = b+n2;

	struct complex sum_a0_a1[KA2_REAL_MAX_AB_SIZE];
	struct complex sum_b0_b1[KA2_REAL_MAX_AB_SIZE];

	sum(a0, a1, sum_a0_a1, n2);
	sum(b0, b1, sum_b0_b1, n2);

	struct complex p1[KA2_REAL_MAX_P_SIZE];
	struct complex p2[KA2_REAL_MAX_P_SIZE];
	struct complex p3[KA2_REAL_MAX_P_SIZE];
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));

	ka2_real_220(a0, b0, p1, n2);
	ka2_real_220(sum_a0_a1, sum_b0_b1, p2, n2);
	ka2_real_220(a1, b1, p3, n2);


	sum_real(c, p1, c, n-1);
	difference_real(c+n2, p3, c+n2, n-1);

	for(int i=n+n2-2; i>=0 ; --i){
		c[i+n2].re -= c[i].re;
	}

	sum_real(c+n2, p2, c+n2, n-1);

}



void b1_real(struct complex a[], struct complex b[], struct complex c[], int n){
	
	#define B1_REAL_MAX_AB_SIZE 218 // n/3
	#define B1_REAL_MAX_P_SIZE 435 // 2*n/3-1

	/* base case */
	if(n == 1){
		multiplication_real(a, b, c, n);
		return;
	}

	/* Assumed n = 3^k for some non-negative integer k */
	int n3 = n/3;
	int p_size = 2*n3-1;


	struct complex *a0 = a;
	struct complex *a1 = a+n3;
	struct complex *a2 = a+2*n3;
	struct complex *b0 = b;
	struct complex *b1 = b+n3;
	struct complex *b2 = b+2*n3;

	struct complex ra0[B1_REAL_MAX_AB_SIZE];
	struct complex ra1[B1_REAL_MAX_AB_SIZE];
	struct complex ra2[B1_REAL_MAX_AB_SIZE];
	struct complex ra3[B1_REAL_MAX_AB_SIZE+2];
	struct complex ra4[B1_REAL_MAX_AB_SIZE+2];
	struct complex rb0[B1_REAL_MAX_AB_SIZE];
	struct complex rb1[B1_REAL_MAX_AB_SIZE];
	struct complex rb2[B1_REAL_MAX_AB_SIZE];
	struct complex rb3[B1_REAL_MAX_AB_SIZE+2];
	struct complex rb4[B1_REAL_MAX_AB_SIZE+2];

	sum_real(a0, a2, ra0, n3);
	sum_real(ra0, a1, ra1, n3);
	difference_real(ra0, a1, ra2, n3);

	ra3[0].re = 0;
	ra3[1].re = a1[0].re;
	for(int i=2; i<n3+1; ++i){
		ra3[i].re = a1[i-1].re + a2[i-2].re;
		
	}
	ra3[n3+1].re = a2[n3-1].re;

	sum_real(a0, ra3, ra4, n3);
	ra4[n3].re = ra3[n3].re;
	ra4[n3+1].re = ra3[n3+1].re;


	sum_real(b0, b2, rb0, n3);
	sum_real(rb0, b1, rb1, n3);
	difference_real(rb0, b1, rb2, n3);

	rb3[0].re =0;
	
	rb3[1].re = b1[0].re;
	for(int i=2; i<n3+1; ++i){
		rb3[i].re = b1[i-1].re + b2[i-2].re;
		
	}
	rb3[n3+1].re = b2[n3-1].re;

	
	sum_real(b0, rb3, rb4, n3);
	rb4[n3].re= rb3[n3].re;
	rb4[n3+1].re = rb3[n3+1].re;


	struct complex p0[B1_REAL_MAX_P_SIZE];
	struct complex p1[B1_REAL_MAX_P_SIZE];
	struct complex p2[B1_REAL_MAX_P_SIZE];
	struct complex p3[B1_REAL_MAX_P_SIZE+4];
	struct complex p4[B1_REAL_MAX_P_SIZE];
	memset(p0, 0, sizeof(p0));
	memset(p1, 0, sizeof(p1));
	memset(p2, 0, sizeof(p2));
	memset(p3, 0, sizeof(p3));
	memset(p4, 0, sizeof(p4));


	ka2_real_218(a0, b0, p0, n3);
	ka2_real_218(ra1, rb1, p1, n3);
	ka2_real_218(ra2, rb2, p2, n3);
	ka2_real_220(ra4, rb4, p3, n3+2);
	ka2_real_218(a2, b2, p4, n3);

	struct complex v0[B1_REAL_MAX_P_SIZE];
	struct complex v1[B1_REAL_MAX_P_SIZE];
	struct complex v2[B1_REAL_MAX_P_SIZE+1];
	struct complex v3[B1_REAL_MAX_P_SIZE+3];
	struct complex v4[B1_REAL_MAX_P_SIZE+2];
	struct complex v[B1_REAL_MAX_P_SIZE+1];
	struct complex u0[B1_REAL_MAX_P_SIZE+1];
	struct complex u[B1_REAL_MAX_P_SIZE+1];

	sum_real(p1, p2, v0, p_size);
	difference_real(p1, p2, v1, p_size);

	v2[0].re= v1[0].re;
	for(int i=1; i<p_size; ++i){
		v2[i].re = v0[i-1].re + v1[i].re;
	}
	v2[p_size].re = v0[p_size-1].re;

	for(int i=0; i<p_size+1; ++i){
		v3[i].re = v2[i].re + p3[i+1].re;
	}
	v3[p_size+1].re = p3[p_size+2].re;
	v3[p_size+2].re = p3[p_size+3].re;

	int sum_re = 0;
	for(int i=p_size+1; i>=0; --i){
		sum_re += v3[i+1].re;

		v4[i].re = sum_re;
	}

	sum_re = 0;
	for(int i=p_size; i>=0; --i){
		sum_re = - sum_re + v4[i+1].re;

		v[i].re = sum_re;
	}


	u0[0].re= v[0].re;
	for(int i=1; i<p_size+1; ++i){
		u0[i].re = v[i].re - p4[i-1].re;
	}


	for(int i=0; i<p_size-1; ++i){
		u[i].re = u0[i].re + p0[i+1].re;
	}
	u[p_size-1].re = u0[p_size-1].re;
	u[p_size].re = u0[p_size].re;


	sum_real(c+4*n3, p4, c+4*n3, p_size);
	sum_real(c+3*n3, u, c+3*n3, p_size+1);
	sum_real(c, p0, c, p_size);
	difference_real(c+n3, u, c+n3, p_size+1);
	difference_real(c+n3, v1, c+n3, p_size);
	difference_real(c+2*n3, p4, c+2*n3, p_size);
	difference_real(c+2*n3, p0, c+2*n3, p_size);
	difference_real(c+2*n3, v0, c+2*n3, p_size);

}


int main(){


	struct complex *a = malloc(654*sizeof(struct complex));  
	struct complex *b = malloc(654*sizeof(struct complex));
	struct complex *c = calloc((2*654-1), sizeof(struct complex));

       int i;
     
//reading a and b from the file inp654test
    FILE *myFile;
    myFile = fopen("inp654test", "r");
    if (myFile == NULL){
        printf("Error Reading File\n");
        exit (0);}

    for (i = 0; i < 654; i++){
        fscanf(myFile, "%d", &a[i].re );
        a[i].im = 0;
    }
    
    for (i = 0; i < 654; i++){
        fscanf(myFile, "%d", &b[i].re );
        b[i].im = 0;
    }
    fclose(myFile);
    

       unsigned long long   LastCycleCount = _rdtsc();
    	clock_t start = clock ();	     
   	for (i=0;i<99999;i++){
           b1_real(a, b, c, 654);
           memset(c, 0, sizeof(struct complex)*(2*654-1));
    	    }
  	b1_real(a, b, c, 654);
	unsigned long long   EndCycleCount = _rdtsc();
        unsigned long long   CyclesElapsed = EndCycleCount - LastCycleCount;
        CyclesElapsed = CyclesElapsed/100000;
       double Timelapsed=(clock()-start)/(double) CLOCKS_PER_SEC;
       Timelapsed=Timelapsed/100000;
       
//writing on c_output file	
	FILE *outFile;
        outFile = fopen("c_output", "w");
        if (outFile == NULL)
	{
    	 printf("Cannot Open File\n");
         exit (0);
         }
      for (i = 0; i < 2*654-1; i++)
       fprintf(outFile, "%d ", ((c[i].re)%3+3)%3 );
       fclose(outFile);
//printing on console
printf("c:");
printf(" ");	
  for( i=0; i<2*654-1; i++)
    printf("%d ", ((c[i].re)%3+3)%3);
    printf("\n\n\n");
    printf("Cycles: %llu\n", CyclesElapsed);
    printf("Multiplication Time:");
    printf( "%f",Timelapsed);
    printf("sec");
    printf("\n \n ");

	free(a);
	free(b);
	free(c);

	return 0;
}
