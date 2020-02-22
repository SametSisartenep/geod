#include <u.h>
#include <libc.h>
#include "geodesic.h"

static uint digits, maxit1, maxit2;
static double ε, realmin, tiny,
		tol0, tol1, tol2, tolb,
		xthresh;

static double
max(double a, double b)
{
	return a > b ? a : b;
}

static double
min(double a, double b)
{
	return a < b ? a : b;
}

static void
swap(double *x, double *y)
{
	double t;

	t = *x;
	*x = *y;
	*y = t;
}

static double
sum(double u, double v, double *t)
{
	double s, up, vpp;

	s = u + v;
	up = s - v;
	vpp = s - up;
	up -= u;
	vpp -= v;
	if(t)
		*t = -(up + vpp);
	/*
	 * error-free sum:
	 * u + v =       s      + t
	 *       = round(u + v) + t
	 */
	return s;
}

static double
log1p(double x)
{
	double y, z;

	y = x + 1;
	z = y - 1;
	/*
	 * Here's the explanation for this magic: y = 1 + z, exactly, and z
	 * approx x, thus log(y)/z (which is nearly constant near z = 0) returns
	 * a good approximation to the true log(1 + x)/x.  The multiplication x *
	 * (log(y)/z) introduces little additional error.
	 */
	return z == 0 ? x : x * log(y)/z;
}

static double
atanh(double x)
{
	double y;

	y = fabs(x);	/* Enforce odd parity */
	y = log1p(2*y/(1 - y))/2;
	return x < 0 ? -y : y;
}

static double
atan2d(double y, double x)
{
	/*
	 * In order to minimize round-off errors, this function rearranges the
	 * arguments so that result of atan2 is in the range [-π/4, π/4] before
	 * converting it to degrees and mapping the result to the correct
	 * quadrant.
	 */
	double ang;
	int q;

	q = 0;
	if(fabs(y) > fabs(x)){
		swap(&x, &y);
		q = 2;
	}
	if(x < 0){
		x = -x;
		++q;
	}
	/* here x >= 0 and x >= abs(y), so angle is in [-pi/4, pi/4] */
	ang = atan2(y, x)/DEG;
	/* Note that atan2d(-0.0, 1.0) will return -0.  However, we expect that
	* atan2d will not be called with y = -0.  If need be, include
	*
	*   case 0: ang = 0 + ang; break;
	*/
	switch(q){
	case 1: ang = (y >= 0 ? 180 : -180) - ang; break;
	case 2: ang =  90 - ang; break;
	case 3: ang = -90 + ang; break;
	}
	return ang;
}

static double
cbrt(double x)
{
	double y;

	y = pow(fabs(x), 1/3.0);
	return x < 0 ? -y : y;
}

static void
norm2(double *sinx, double *cosx)
{
	double r;

	r = hypot(*sinx, *cosx);
	*sinx /= r;
	*cosx /= r;
}

static double
copysign(double x, double y)
{
	return fabs(x)*(y < 0 || (y == 0 && 1/y < 0) ? -1 : 1);
}

static double
AngNormalize(double x)
{
	double y;

	y = fmod(x, 360.0);
	return y <= -180 ? y + 360 :
		y <= 180 ? y : y - 360;
}

static double
AngRound(double x)
{
	double y, z;

	z = 1/16.0;
	if(x == 0)
		return 0;
	y = fabs(x);
	/* The compiler mustn't "simplify" z - (z - y) to y */
	y = y < z ? z - (z - y) : y;
	return x < 0 ? -y : y;
}

static double
AngDiff(double x, double y, double *e)
{
	double t, d;

	d = AngNormalize(sum(AngNormalize(-x), AngNormalize(y), &t));
	/*
	 * Here y - x = d + t (mod 360), exactly, where d is in (-180,180] and
	 * abs(t) <= eps (eps = 2^-45 for doubles).  The only case where the
	 * addition of t takes the result outside the range (-180,180] is d = 180
	 * and t > 0.  The case, d = -180 + eps, t = -eps, can't happen, since
	 * sum would have returned the exact result in such a case (i.e., given t
	 * = 0).
	 */
	return sum(d == 180 && t > 0 ? -180 : d, t, e);
}

static double
LatFix(double x)
{
	return fabs(x) > 90 ? NaN() : x;
}

static void
sincosd(double x, double *sinx, double *cosx)
{
	/*
	 * In order to minimize round-off errors, this function exactly reduces
	 * the argument to the range [-45, 45] before converting it to radians.
	 */
	double r, s, c;
	int q;

	r = fmod(x, 360.0);
	q = !isNaN(r) ? (int)floor(r/90 + 0.5) : 0;
	r -= 90 * q;
	/* now |r| <= 45 */
	s = sin(r*DEG);
	c = cos(r*DEG);
	switch(q & 3){
	case 0: *sinx =  s; *cosx =  c; break;
	case 1: *sinx =  c; *cosx = -s; break;
	case 2: *sinx = -s; *cosx = -c; break;
	case 3: *sinx = -c; *cosx =  s; break;
	}
	if(x != 0){
		*sinx += 0.0;
		*cosx += 0.0;
	}
}

static double
polyval(int N, double p[], double x)
{
	double y;

	y = N < 0 ? 0 : *p++;
	while(--N >= 0)
		y *= x + *p++;
	return y;
}

/* The scale factor A3 = mean value of (d/dsigma)I3 */
static void
A3coeff(Geodesic *g)
{
	static double coeff[] = {
		/* A3, coeff of eps^5, polynomial in n of order 0 */
		-3, 128,
		/* A3, coeff of eps^4, polynomial in n of order 1 */
		-2, -3, 64,
		/* A3, coeff of eps^3, polynomial in n of order 2 */
		-1, -3, -1, 16,
		/* A3, coeff of eps^2, polynomial in n of order 2 */
		3, -1, -2, 8,
		/* A3, coeff of eps^1, polynomial in n of order 1 */
		1, -1, 2,
		/* A3, coeff of eps^0, polynomial in n of order 0 */
		1, 1,
	};
	int o, k, j, m;

	o = k = 0;
	for(j = nA3 - 1; j >= 0; --j){	/* coeff of eps^j */
		m = min(nA3 - j-1, j);	/* order of polynomial in n */
		g->A3x[k++] = polyval(m, coeff+o, g->n) / coeff[o+m+1];
		o += m + 2;
	}
}

/* The coefficients C3[l] in the Fourier expansion of B3 */
static void
C3coeff(Geodesic *g)
{
	static double coeff[] = {
		/* C3[1], coeff of eps^5, polynomial in n of order 0 */
		3, 128,
		/* C3[1], coeff of eps^4, polynomial in n of order 1 */
		2, 5, 128,
		/* C3[1], coeff of eps^3, polynomial in n of order 2 */
		-1, 3, 3, 64,
		/* C3[1], coeff of eps^2, polynomial in n of order 2 */
		-1, 0, 1, 8,
		/* C3[1], coeff of eps^1, polynomial in n of order 1 */
		-1, 1, 4,
		/* C3[2], coeff of eps^5, polynomial in n of order 0 */
		5, 256,
		/* C3[2], coeff of eps^4, polynomial in n of order 1 */
		1, 3, 128,
		/* C3[2], coeff of eps^3, polynomial in n of order 2 */
		-3, -2, 3, 64,
		/* C3[2], coeff of eps^2, polynomial in n of order 2 */
		1, -3, 2, 32,
		/* C3[3], coeff of eps^5, polynomial in n of order 0 */
		7, 512,
		/* C3[3], coeff of eps^4, polynomial in n of order 1 */
		-10, 9, 384,
		/* C3[3], coeff of eps^3, polynomial in n of order 2 */
		5, -9, 5, 192,
		/* C3[4], coeff of eps^5, polynomial in n of order 0 */
		7, 512,
		/* C3[4], coeff of eps^4, polynomial in n of order 1 */
		-14, 7, 512,
		/* C3[5], coeff of eps^5, polynomial in n of order 0 */
		21, 2560,
	};
	int o, k, l, j, m;

	o = k = 0;
	for(l = 1; l < nC3; ++l){	/* l is index of C3[l] */
		for(j = nC3 - 1; j >= l; --j){	/* coeff of eps^j */
			m = min(nC3 - j-1, j); /* order of polynomial in n */
			g->C3x[k++] = polyval(m, coeff+o, g->n) / coeff[o+m+1];
			o += m + 2;
		}
	}
}

/* The coefficients C4[l] in the Fourier expansion of I4 */
static void
C4coeff(Geodesic *g)
{
	static double coeff[] = {
		/* C4[0], coeff of eps^5, polynomial in n of order 0 */
		97, 15015,
		/* C4[0], coeff of eps^4, polynomial in n of order 1 */
		1088, 156, 45045,
		/* C4[0], coeff of eps^3, polynomial in n of order 2 */
		-224, -4784, 1573, 45045,
		/* C4[0], coeff of eps^2, polynomial in n of order 3 */
		-10656, 14144, -4576, -858, 45045,
		/* C4[0], coeff of eps^1, polynomial in n of order 4 */
		64, 624, -4576, 6864, -3003, 15015,
		/* C4[0], coeff of eps^0, polynomial in n of order 5 */
		100, 208, 572, 3432, -12012, 30030, 45045,
		/* C4[1], coeff of eps^5, polynomial in n of order 0 */
		1, 9009,
		/* C4[1], coeff of eps^4, polynomial in n of order 1 */
		-2944, 468, 135135,
		/* C4[1], coeff of eps^3, polynomial in n of order 2 */
		5792, 1040, -1287, 135135,
		/* C4[1], coeff of eps^2, polynomial in n of order 3 */
		5952, -11648, 9152, -2574, 135135,
		/* C4[1], coeff of eps^1, polynomial in n of order 4 */
		-64, -624, 4576, -6864, 3003, 135135,
		/* C4[2], coeff of eps^5, polynomial in n of order 0 */
		8, 10725,
		/* C4[2], coeff of eps^4, polynomial in n of order 1 */
		1856, -936, 225225,
		/* C4[2], coeff of eps^3, polynomial in n of order 2 */
		-8448, 4992, -1144, 225225,
		/* C4[2], coeff of eps^2, polynomial in n of order 3 */
		-1440, 4160, -4576, 1716, 225225,
		/* C4[3], coeff of eps^5, polynomial in n of order 0 */
		-136, 63063,
		/* C4[3], coeff of eps^4, polynomial in n of order 1 */
		1024, -208, 105105,
		/* C4[3], coeff of eps^3, polynomial in n of order 2 */
		3584, -3328, 1144, 315315,
		/* C4[4], coeff of eps^5, polynomial in n of order 0 */
		-128, 135135,
		/* C4[4], coeff of eps^4, polynomial in n of order 1 */
		-2560, 832, 405405,
		/* C4[5], coeff of eps^5, polynomial in n of order 0 */
		128, 99099,
	};
	int o, k, l, j, m;

	o = k = 0;
	for(l = 0; l < nC4; ++l){	/* l is index of C4[l] */
		for(j = nC4 - 1; j >= l; --j){	/* coeff of eps^j */
			m = nC4 - j-1;	/* order of polynomial in n */
			g->C4x[k++] = polyval(m, coeff+o, g->n) / coeff[o+m+1];
			o += m + 2;
		}
	}
}

/* The scale factor A1-1 = mean value of (d/dsigma)I1 - 1 */
static double
A1m1f(double eps)
{
	static double coeff[] = {
		/* (1-eps)*A1-1, polynomial in eps2 of order 3 */
		1, 4, 64, 0, 256,
	};
	double t;
	int m;

	m = nA1/2;
	t = polyval(m, coeff, eps*eps)/coeff[m+1];
	return (t + eps)/(1 - eps);
}

/* The coefficients C1[l] in the Fourier expansion of B1 */
static void
C1f(double eps, double c[])
{
	static double coeff[] = {
		/* C1[1]/eps^1, polynomial in eps2 of order 2 */
		-1, 6, -16, 32,
		/* C1[2]/eps^2, polynomial in eps2 of order 2 */
		-9, 64, -128, 2048,
		/* C1[3]/eps^3, polynomial in eps2 of order 1 */
		9, -16, 768,
		/* C1[4]/eps^4, polynomial in eps2 of order 1 */
		3, -5, 512,
		/* C1[5]/eps^5, polynomial in eps2 of order 0 */
		-7, 1280,
		/* C1[6]/eps^6, polynomial in eps2 of order 0 */
		-7, 2048,
	};
	double eps², d;
	int o, l, m;

	eps² = eps*eps;
	d = eps;
	o = 0;
	for(l = 1; l <= nC1; ++l){	/* l is index of C1p[l] */
		m = (nC1 - l)/2;	/* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeff+o, eps²)/coeff[o+m+1];
		o += m + 2;
		d *= eps;
	}
}

/* The coefficients C1p[l] in the Fourier expansion of B1p */
static void
C1pf(double eps, double c[])
{
	static double coeff[] = {
		/* C1p[1]/eps^1, polynomial in eps2 of order 2 */
		205, -432, 768, 1536,
		/* C1p[2]/eps^2, polynomial in eps2 of order 2 */
		4005, -4736, 3840, 12288,
		/* C1p[3]/eps^3, polynomial in eps2 of order 1 */
		-225, 116, 384,
		/* C1p[4]/eps^4, polynomial in eps2 of order 1 */
		-7173, 2695, 7680,
		/* C1p[5]/eps^5, polynomial in eps2 of order 0 */
		3467, 7680,
		/* C1p[6]/eps^6, polynomial in eps2 of order 0 */
		38081, 61440,
	};
	double eps², d;
	int o, l, m;
	
	eps² = eps*eps,
	d = eps;
	o = 0;
	for(l = 1; l <= nC1p; ++l){	/* l is index of C1p[l] */
		m = (nC1p - l)/2;	/* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeff+o, eps²)/coeff[o+m+1];
		o += m + 2;
		d *= eps;
	}
}

/* The scale factor A2-1 = mean value of (d/dsigma)I2 - 1 */
static double
A2m1f(double eps)
{
	static double coeff[] = {
		/* (eps+1)*A2-1, polynomial in eps2 of order 3 */
		-11, -28, -192, 0, 256,
	};
	double t;
	int m;

	m = nA2/2;
	t = polyval(m, coeff, eps*eps)/coeff[m+1];
	return (t - eps)/(1 + eps);
}

/* The coefficients C2[l] in the Fourier expansion of B2 */
static void
C2f(double eps, double c[])
{
	static double coeff[] = {
		/* C2[1]/eps^1, polynomial in eps2 of order 2 */
		1, 2, 16, 32,
		/* C2[2]/eps^2, polynomial in eps2 of order 2 */
		35, 64, 384, 2048,
		/* C2[3]/eps^3, polynomial in eps2 of order 1 */
		15, 80, 768,
		/* C2[4]/eps^4, polynomial in eps2 of order 1 */
		7, 35, 512,
		/* C2[5]/eps^5, polynomial in eps2 of order 0 */
		63, 1280,
		/* C2[6]/eps^6, polynomial in eps2 of order 0 */
		77, 2048,
	};
	double eps², d;
	int o, l, m;
	
	eps² = eps*eps,
	d = eps;
	o = 0;
	for(l = 1; l <= nC2; ++l){	/* l is index of C2[l] */
		m = (nC2 - l)/2;	/* order of polynomial in eps^2 */
		c[l] = d * polyval(m, coeff+o, eps²)/coeff[o+m+1];
		o += m + 2;
		d *= eps;
	}
}

static double
A3f(Geodesic *g, double eps)
{
	/* Evaluate A3 */
	return polyval(nA3 - 1, g->A3x, eps);
}

static void
C3f(Geodesic *g, double eps, double c[])
{
	/*
	 * Evaluate C3 coeffs
	 * Elements c[1] through c[nC3 - 1] are set
	 */
	double mult;
	int o, l, m;

	mult = 1;
	o = 0;
	for(l = 1; l < nC3; ++l){	/* l is index of C3[l] */
		m = nC3 - l-1;	/* order of polynomial in eps */
		mult *= eps;
		c[l] = mult*polyval(m, g->C3x + o, eps);
		o += m + 1;
	}
}

static void
C4f(Geodesic *g, double eps, double c[])
{
	/*
	 * Evaluate C4 coeffs
	 * Elements c[0] through c[nC4 - 1] are set
	 */
	double mult;
	int o, l, m;

	mult = 1;
	o = 0;
	for(l = 0; l < nC4; ++l){	/* l is index of C4[l] */
		m = nC4 - l-1;	/* order of polynomial in eps */
		c[l] = mult*polyval(m, g->C4x + o, eps);
		o += m + 1;
		mult *= eps;
	}
}

static double
SinCosSeries(int sinp, double sinx, double cosx, double c[], int n)
{
	/*
	 * Evaluate
	 * y = sinp ? sum(c[i] * sin( 2*i    * x), i, 1, n) :
	 *            sum(c[i] * cos((2*i+1) * x), i, 0, n-1)
	 * using Clenshaw summation.  N.B. c[0] is unused for sin series
	 * Approx operation count = (n + 5) mult and (2*n + 2) add
	 */
	double ar, y0, y1;

	c += (n + sinp);	/* Point to one beyond last element */
	ar = 2*(cosx - sinx)*(cosx + sinx);	/* 2 * cos(2 * x) */
	/* accumulators for sum */
	y0 = (n & 1) ? *--c : 0;
	y1 = 0;
	/* Now n is even */
	n /= 2;
	while(n--){
		/* Unroll loop x 2, so accumulators return to their original role */
		y1 = ar * y0 - y1 + *--c;
		y0 = ar * y1 - y0 + *--c;
	}
	return sinp ?
		2 * sinx * cosx * y0 :      /* sin(2 * x) * y0 */
		cosx * (y0 - y1);         /* cos(x) * (y0 - y1) */
}

static void
Lengths(Geodesic *g,
	double eps, double sig12,
	double ssig1, double csig1, double dn1,
	double ssig2, double csig2, double dn2,
	double cbet1, double cbet2,
	double *ps12b, double *pm12b, double *pm0,
	double *pM12, double *pM21,
	/* Scratch area of the right size */
	double Ca[])
{
	double m0, J12, A1, A2;
	double Cb[nC];
	int redlp;

	m0 = J12 = A1 = A2 = 0;

	/* Return m12b = (reduced length)/b; also calculate s12b = distance/b,
	* and m0 = coefficient of secular term in expression for reduced length. */
	redlp = pm12b || pm0 || pM12 || pM21;
	if(ps12b || redlp){
		A1 = A1m1f(eps);
		C1f(eps, Ca);
		if(redlp){
			A2 = A2m1f(eps);
			C2f(eps, Cb);
			m0 = A1 - A2;
			A2 = 1 + A2;
		}
		A1 = 1 + A1;
	}
	if(ps12b){
		double B1;

		B1 = SinCosSeries(1, ssig2, csig2, Ca, nC1) -
			SinCosSeries(1, ssig1, csig1, Ca, nC1);
		/* Missing a factor of b */
		*ps12b = A1*(sig12 + B1);
		if(redlp){
			double B2;

			B2 = SinCosSeries(1, ssig2, csig2, Cb, nC2) -
				SinCosSeries(1, ssig1, csig1, Cb, nC2);
			J12 = m0 * sig12 + (A1*B1 - A2*B2);
		}
	}else if(redlp){
		/* Assume here that nC1 >= nC2 */
		int l;

		for(l = 1; l <= nC2; ++l)
			Cb[l] = A1*Ca[l] - A2*Cb[l];
		J12 = m0 * sig12 + (SinCosSeries(1, ssig2, csig2, Cb, nC2) -
			SinCosSeries(1, ssig1, csig1, Cb, nC2));
	}
	if(pm0)
		*pm0 = m0;
	if(pm12b)
		/*
		 * Missing a factor of b.
		 * Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
		 * accurate cancellation in the case of coincident points.
		 */
		*pm12b = dn2 * (csig1*ssig2) - dn1 * (ssig1*csig2) - csig1*csig2 * J12;
	if(pM12 || pM21){
		double csig12, t;

		csig12 = csig1*csig2 + ssig1*ssig2;
		t = g->ep2*(cbet1 - cbet2)*(cbet1 + cbet2)/(dn1 + dn2);
		if(pM12)
			*pM12 = csig12 + (t * ssig2 - csig2 * J12) * ssig1 / dn1;
		if(pM21)
			*pM21 = csig12 - (t * ssig1 - csig1 * J12) * ssig2 / dn2;
	}
}

static double
Astroid(double x, double y)
{
	/* Solve k⁴+2*k³-(x²+y²-1)*k²-2*y²*k-y² = 0 for positive root k. */
	double k, p, q, r;

	p = x*x;
	q = y*y;
	r = (p + q - 1)/6;
	if(!(q == 0 && r <= 0)){
		double S, r2, r3, disc;
		double u, v, uv, w;

		/*
		 * Avoid possible division by zero when r = 0 by multiplying equations
		 * for s and t by r^3 and r, resp.
		 */
		S = p*q / 4;	/* S = r^3 * s */
		r2 = r*r;
		r3 = r*r2;
		/* The discriminant of the quadratic equation for T3.  This is zero on
		* the evolute curve p^(1/3)+q^(1/3) = 1 */
		disc = S*(S + 2*r3);
		u = r;
		if(disc >= 0){
			double T3, T;

			T3 = S + r3;
			/*
			 * Pick the sign on the sqrt to maximize abs(T3).  This minimizes loss
			 * of precision due to cancellation.  The result is unchanged because
			 * of the way the T is used in definition of u.
			 */
			T3 += T3 < 0 ? -sqrt(disc) : sqrt(disc);	/* T3 = (r * t)^3 */
			/* N.B. cbrtx always returns the real root.  cbrtx(-8) = -2. */
			T = cbrt(T3);	/* T = r * t */
			/* T can be zero; but then r2 / T -> 0. */
			u += T + (T != 0 ? r2/T : 0);
		}else{
			double ang;

			/* T is complex, but the way u is defined the result is real. */
			ang = atan2(sqrt(-disc), -(S + r3));
			/*
			 * There are three possible cube roots.  We choose the root which
			 * avoids cancellation.  Note that disc < 0 implies that r < 0.
			 */
			u += 2*r*cos(ang/3);
		}
		v = sqrt(u*u + q);	/* guaranteed positive */
		/* Avoid loss of accuracy when u < 0. */
		uv = u < 0 ? q/(v - u) : u + v;	/* u+v, guaranteed positive */

		w = (uv - q)/(2*v);	/* positive? */
		/*
		 * Rearrange expression for k to avoid loss of accuracy due to
		 * subtraction.  Division by 0 not possible because uv > 0, w >= 0.
		 */
		k = uv/(sqrt(uv + w*w) + w);	/* guaranteed positive */
	}else{	/* q == 0 && r <= 0 */
		/*
		 * y = 0 with |x| <= 1.  Handle this case directly.
		 * for y small, positive root is k = abs(y)/sqrt(1-x^2)
		 */
		k = 0;
	}
	return k;
}

static double
InverseStart(Geodesic *g,
	double sbet1, double cbet1, double dn1,
	double sbet2, double cbet2, double dn2,
	double lam12, double slam12, double clam12,
	double *psalp1, double *pcalp1,
	/* Only updated if return val >= 0 */
	double *psalp2, double *pcalp2,
	/* Only updated for short lines */
	double *pdnm,
	/* Scratch area of the right size */
	double Ca[])
{
	double salp1, calp1, salp2, calp2, dnm;
	double sig12, sbet12, cbet12, sbet12a;
	double somg12, comg12, ssig12, csig12;
	int shortline;

	salp1 = calp1 = salp2 = calp2 = dnm = 0;
	/*
	 * Return a starting point for Newton's method in salp1 and calp1 (function
	 * value is -1).  If Newton's method doesn't need to be used, return also
	 * salp2 and calp2 and function value is sig12.
	 */
	sig12 = -1;	/* Return value */
	/* bet12 = bet2 - bet1 in [0, π); bet12a = bet2 + bet1 in (-π, 0] */
	sbet12 = sbet2*cbet1 - cbet2*sbet1;
	cbet12 = cbet2*cbet1 + sbet2*sbet1;
	shortline = cbet12 >= 0 && sbet12 < 0.5 &&
		cbet2*lam12 < 0.5;
	sbet12a = sbet2 * cbet1 + cbet2 * sbet1;
	if(shortline){
		double sbetm2, omg12;

		sbetm2 = (sbet1 + sbet2)*(sbet1 + sbet2);
		/*
		 * sin((bet1+bet2)/2)^2
		 * =  (sbet1 + sbet2)^2 / ((sbet1 + sbet2)^2 + (cbet1 + cbet2)^2)
		 */
		sbetm2 /= sbetm2 + (cbet1 + cbet2)*(cbet1 + cbet2);
		dnm = sqrt(1 + g->ep2 * sbetm2);
		omg12 = lam12 / (g->f1 * dnm);
		somg12 = sin(omg12);
		comg12 = cos(omg12);
	}else{
		somg12 = slam12;
		comg12 = clam12;
	}

	salp1 = cbet2*somg12;
	calp1 = comg12 >= 0 ?
		sbet12 + cbet2 * sbet1 * somg12*somg12 / (1 + comg12) :
		sbet12a - cbet2 * sbet1 * somg12*somg12 / (1 - comg12);

	ssig12 = hypot(salp1, calp1);
	csig12 = sbet1 * sbet2 + cbet1 * cbet2 * comg12;

	if(shortline && ssig12 < g->etol2){
		/* really short lines */
		salp2 = cbet1 * somg12;
		calp2 = sbet12 - cbet1 * sbet2 *
			(comg12 >= 0 ? somg12*somg12 / (1 + comg12) : 1 - comg12);
		norm2(&salp2, &calp2);
		/* Set return value */
		sig12 = atan2(ssig12, csig12);
	}else if(fabs(g->n) > 0.1 || /* No astroid calc if too eccentric */
		csig12 >= 0 ||
		ssig12 >= 6 * fabs(g->n) * PI * cbet1*cbet1){
		/* Nothing to do, zeroth order spherical approximation is OK */
	}else{
		/*
		 * Scale lam12 and bet2 to x, y coordinate system where antipodal point
		 * is at origin and singular point is at y = 0, x = -1.
		 */
		double x, y, lamscale, betscale;
		double lam12x;

		lam12x = atan2(-slam12, -clam12); /* lam12 - π */
		if(g->f >= 0){	/* In fact f == 0 does not get here */
			/* x = dlong, y = dlat */
			{
				double k2, eps;
	
				k2 = sbet1*sbet1 * g->ep2;
				eps = k2 / (2 * (1 + sqrt(1 + k2)) + k2);
				lamscale = g->f * cbet1 * A3f(g, eps) * PI;
			}
			betscale = lamscale * cbet1;
			x = lam12x / lamscale;
			y = sbet12a / betscale;
		}else{	/* f < 0 */
	
			/* x = dlat, y = dlong */
			double cbet12a, bet12a;
			double m12b, m0;
	
			cbet12a = cbet2*cbet1 - sbet2*sbet1;
			bet12a = atan2(sbet12a, cbet12a);
			/*
			 * In the case of lon12 = 180, this repeats a calculation made in
			 * Inverse.
			 */
			Lengths(g, g->n, PI + bet12a,
				sbet1, -cbet1, dn1, sbet2, cbet2, dn2,
				cbet1, cbet2, nil, &m12b, &m0, nil, nil, Ca);
			x = -1 + m12b / (cbet1*cbet2*m0 * PI);
			betscale = x < -0.01 ? sbet12a/x :
				-g->f * cbet1*cbet1 * PI;
			lamscale = betscale / cbet1;
			y = lam12x / lamscale;
		}

		if(y > -tol1 && x > -1 - xthresh){
			/* strip near cut */
			if(g->f >= 0){
				salp1 = min(1, -x);
				calp1 = -sqrt(1 - salp1*salp1);
			}else{
				calp1 = max(x > -tol1 ? 0 : -1, x);
				salp1 = sqrt(1 - calp1*calp1);
			}
		}else{
			/*
			 * Estimate alp1, by solving the astroid problem.
			 *
			 * Could estimate alpha1 = theta + pi/2, directly, i.e.,
			 *   calp1 = y/k; salp1 = -x/(1+k);  for f >= 0
			 *   calp1 = x/(1+k); salp1 = -y/k;  for f < 0 (need to check)
			 *
			 * However, it's better to estimate omg12 from astroid and use
			 * spherical formula to compute alp1.  This reduces the mean number of
			 * Newton iterations for astroid cases from 2.24 (min 0, max 6) to 2.12
			 * (min 0 max 5).  The changes in the number of iterations are as
			 * follows:
			 *
			 * change percent
			 *    1       5
			 *    0      78
			 *   -1      16
			 *   -2       0.6
			 *   -3       0.04
			 *   -4       0.002
			 *
			 * The histogram of iterations is (m = number of iterations estimating
			 * alp1 directly, n = number of iterations estimating via omg12, total
			 * number of trials = 148605):
			 *
			 *  iter    m      n
			 *    0   148    186
			 *    1 13046  13845
			 *    2 93315 102225
	
			 *    3 36189  32341
			 *    4  5396      7
			 *    5   455      1
			 *    6    56      0
			 *
			 * Because omg12 is near π, estimate work with omg12a = π - omg12
			 */
			double k, omg12a;

			k = Astroid(x, y);
			omg12a = lamscale*(g->f >= 0 ? -x * k/(1 + k) : -y * (1 + k)/k);
			somg12 = sin(omg12a);
			comg12 = -cos(omg12a);
			/* Update spherical estimate of alp1 using omg12 instead of lam12 */
			salp1 = cbet2*somg12;
			calp1 = sbet12a - cbet2 * sbet1 * somg12*somg12/(1 - comg12);
		}
	}
	/* Sanity check on starting guess.  Backwards check allows NaN through. */
	if(!(salp1 <= 0))
		norm2(&salp1, &calp1);
	else{
		salp1 = 1;
		calp1 = 0;
	}

	*psalp1 = salp1;
	*pcalp1 = calp1;
	if(shortline)
		*pdnm = dnm;
	if(sig12 >= 0){
		*psalp2 = salp2;
		*pcalp2 = calp2;
	}
	return sig12;
}

static double
Lambda12(Geodesic *g,
	double sbet1, double cbet1, double dn1,
	double sbet2, double cbet2, double dn2,
	double salp1, double calp1,
	double slam120, double clam120,
	double *psalp2, double *pcalp2,
	double *psig12,
	double *pssig1, double *pcsig1,
	double *pssig2, double *pcsig2,
	double *peps,
	double *pdomg12,
	int diffp, double *pdlam12,
	/* Scratch area of the right size */
	double Ca[])
{
	double salp2, calp2, sig12, ssig1, csig1, ssig2, csig2, eps, domg12, dlam12;
	double salp0, calp0;
	double somg1, comg1, somg2, comg2, somg12, comg12, lam12;
	double B312, eta, k2;

	ssig1 = csig1 = ssig2 = csig2 = dlam12 = 0;
	if(sbet1 == 0 && calp1 == 0)
		/*
		 * Break degeneracy of equatorial line.  This case has already been
		 * handled.
		 */
		calp1 = -tiny;

	/* sin(alp1) * cos(bet1) = sin(alp0) */
	salp0 = salp1 * cbet1;
	calp0 = hypot(calp1, salp1*sbet1); /* calp0 > 0 */

	/*
	 * tan(bet1) = tan(sig1) * cos(alp1)
	 * tan(omg1) = sin(alp0) * tan(sig1) = tan(omg1)=tan(alp1)*sin(bet1)
	 */
	ssig1 = sbet1;
	somg1 = salp0 * sbet1;
	csig1 = comg1 = calp1 * cbet1;
	norm2(&ssig1, &csig1);
	/* norm2(&somg1, &comg1); -- don't need to normalize! */

	/* Enforce symmetries in the case abs(bet2) = -bet1.  Need to be careful
	* about this case, since this can yield singularities in the Newton
	* iteration.
	* sin(alp2) * cos(bet2) = sin(alp0) */
	salp2 = cbet2 != cbet1 ? salp0/cbet2 : salp1;
	/*
	 * calp2 = sqrt(1 - sq(salp2))
	 *       = sqrt(sq(calp0) - sq(sbet2)) / cbet2
	 * and subst for calp0 and rearrange to give (choose positive sqrt
	 * to give alp2 in [0, pi/2]).
	 */
	calp2 = cbet2 != cbet1 || fabs(sbet2) != -sbet1 ?
		sqrt((calp1 * cbet1)*(calp1 * cbet1) +
		(cbet1 < -sbet1 ? (cbet2 - cbet1) * (cbet1 + cbet2) :
		(sbet1 - sbet2) * (sbet1 + sbet2))) / cbet2 :
		fabs(calp1);
	/*
	 * tan(bet2) = tan(sig2) * cos(alp2)
	 * tan(omg2) = sin(alp0) * tan(sig2).
	 */
	ssig2 = sbet2;
	somg2 = salp0 * sbet2;
	csig2 = comg2 = calp2 * cbet2;
	norm2(&ssig2, &csig2);
	/* norm2(&somg2, &comg2); -- don't need to normalize! */

	/* sig12 = sig2 - sig1, limit to [0, pi] */
	sig12 = atan2(max(0, csig1*ssig2 - ssig1*csig2), csig1*csig2 + ssig1*ssig2);

	/* omg12 = omg2 - omg1, limit to [0, pi] */
	somg12 = max(0, comg1*somg2 - somg1*comg2);
	comg12 = comg1*comg2 + somg1*somg2;
	/* eta = omg12 - lam120 */
	eta = atan2(somg12*clam120 - comg12*slam120, comg12*clam120 + somg12*slam120);
	k2 = calp0*calp0 * g->ep2;
	eps = k2/(2*(1 + sqrt(1 + k2)) + k2);
	C3f(g, eps, Ca);
	B312 = SinCosSeries(1, ssig2, csig2, Ca, nC3-1) -
		SinCosSeries(1, ssig1, csig1, Ca, nC3-1);
	domg12 = -g->f * A3f(g, eps) * salp0 * (sig12 + B312);
	lam12 = eta + domg12;

	if(diffp){
		if(calp2 == 0)
			dlam12 = -2*g->f1*dn1 / sbet1;
		else{
			Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
				cbet1, cbet2, nil, &dlam12, nil, nil, nil, Ca);
			dlam12 *= g->f1/(calp2*cbet2);
		}
	}

	*psalp2 = salp2;
	*pcalp2 = calp2;
	*psig12 = sig12;
	*pssig1 = ssig1;
	*pcsig1 = csig1;
	*pssig2 = ssig2;
	*pcsig2 = csig2;
	*peps = eps;
	*pdomg12 = domg12;
	if(diffp)
		*pdlam12 = dlam12;
	return lam12;
}


static void
init(void)
{
	digits = 53;
	ε = pow(0.5, digits - 1);
	realmin = pow(0.5, 1022);
	maxit1 = 20;
	maxit2 = maxit1 + digits + 10;
	tiny = sqrt(realmin);
	tol0 = ε;
	/*
	 * Increase multiplier in defn of tol1 from 100 to 200 to fix inverse case
	 * 52.784459512564 0 -52.784459512563990912 179.634407464943777557
	 * which otherwise failed for Visual Studio 10 (Release and Debug)
	 */
	tol1 = 200 * tol0;
	tol2 = sqrt(tol0);
	/* Check on bisection interval */
	tolb = tol0 * tol2;
	xthresh = 1000 * tol2;
}

static void
initgeodline_int(Geodline *l, Geodesic *g,
	double lat1, double lon1,
	double azi1, double salp1, double calp1,
	uint caps)
{
	double cbet1, sbet1, eps;

	l->a = g->a;
	l->f = g->f;
	l->b = g->b;
	l->c2 = g->c2;
	l->f1 = g->f1;
	/* If caps is 0 assume the standard direct calculation */
	l->caps = (caps ? caps : GDistanceIn | GLongitude) |
	/* always allow latitude and azimuth and unrolling of longitude */
		GLatitude | GAzimuth | GUnrollLon;
	l->lat1 = LatFix(lat1);
	l->lon1 = lon1;
	l->azi1 = azi1;
	l->salp1 = salp1;
	l->calp1 = calp1;

	sincosd(AngRound(l->lat1), &sbet1, &cbet1);
	sbet1 *= l->f1;
	/* Ensure cbet1 = +epsilon at poles */
	norm2(&sbet1, &cbet1);
	cbet1 = max(tiny, cbet1);
	l->dn1 = sqrt(1 + g->ep2 * sbet1*sbet1);

	/* Evaluate alp0 from sin(alp1) * cos(bet1) = sin(alp0), */
	l->salp0 = l->salp1 * cbet1; /* alp0 in [0, π/2 - |bet1|] */
	/*
	 * Alt: calp0 = hypot(sbet1, calp1 * cbet1).  The following
	 * is slightly better (consider the case salp1 = 0).
	 */
	l->calp0 = hypot(l->calp1, l->salp1 * sbet1);
	/*
	 * Evaluate sig with tan(bet1) = tan(sig1) * cos(alp1).
	 * sig = 0 is nearest northward crossing of equator.
	 * With bet1 = 0, alp1 = pi/2, we have sig1 = 0 (equatorial line).
	 * With bet1 =  pi/2, alp1 = -pi, sig1 =  pi/2
	 * With bet1 = -pi/2, alp1 =  0 , sig1 = -pi/2
	 * Evaluate omg1 with tan(omg1) = sin(alp0) * tan(sig1).
	 * With alp0 in (0, pi/2], quadrants for sig and omg coincide.
	 * No atan2(0,0) ambiguity at poles since cbet1 = +epsilon.
	 * With alp0 = 0, omg1 = 0 for alp1 = 0, omg1 = pi for alp1 = pi.
	 */
	l->ssig1 = sbet1;
	l->somg1 = l->salp0 * sbet1;
	l->csig1 = l->comg1 = sbet1 != 0 || l->calp1 != 0 ? cbet1 * l->calp1 : 1;
	norm2(&l->ssig1, &l->csig1); /* sig1 in (-π, π] */
	/* norm2(somg1, comg1); -- don't need to normalize! */

	l->k2 = l->calp0*l->calp0 * g->ep2;
	eps = l->k2 / (2*(1 + sqrt(1 + l->k2)) + l->k2);

	if(l->caps & CapC1){
		double s, c;

		l->A1m1 = A1m1f(eps);
		C1f(eps, l->C1a);
		l->B11 = SinCosSeries(1, l->ssig1, l->csig1, l->C1a, nC1);
		s = sin(l->B11);
		c = cos(l->B11);
		/* tau1 = sig1 + B11 */
		l->stau1 = l->ssig1*c + l->csig1*s;
		l->ctau1 = l->csig1*c - l->ssig1*s;
		/*
		 * Not necessary because C1pa reverts C1a
		 * B11 = -SinCosSeries(1, stau1, ctau1, C1pa, nC1p);
		 */
	}
	if(l->caps & CapC1p)
		C1pf(eps, l->C1pa);
	if(l->caps & CapC2){
		l->A2m1 = A2m1f(eps);
		C2f(eps, l->C2a);
		l->B21 = SinCosSeries(1, l->ssig1, l->csig1, l->C2a, nC2);
	}
	if(l->caps & CapC3){
		C3f(g, eps, l->C3a);
		l->A3c = -l->f * l->salp0 * A3f(g, eps);
		l->B31 = SinCosSeries(1, l->ssig1, l->csig1, l->C3a, nC3-1);
	}
	if(l->caps & CapC4){
		C4f(g, eps, l->C4a);
		/* Multiplier = a^2 * e^2 * cos(alpha0) * sin(alpha0) */
		l->A4 = l->a*l->a * l->calp0 * l->salp0 * g->e2;
		l->B41 = SinCosSeries(0, l->ssig1, l->csig1, l->C4a, nC4);
	}

	l->a13 = l->s13 = NaN();
}

void
initgeodline(Geodline *l, Geodesic *g, double lat1, double lon1, double azi1, uint caps)
{
	double salp1, calp1;

	azi1 = AngNormalize(azi1);
	/* Guard against underflow in salp0 */
	sincosd(AngRound(azi1), &salp1, &calp1);
	initgeodline_int(l, g, lat1, lon1, azi1, salp1, calp1, caps);
}

double
geod_genposition(Geodline *l,
	uint flags, double s12_a12,
	double *plat2, double *plon2, double *pazi2,
	double *ps12, double *pm12,
	double *pM12, double *pM21,
	double *pS12)
{
	double lat2, lon2, azi2, s12, m12, M12, M21, S12;
	double sig12, ssig12, csig12, B12, AB1;
	double omg12, lam12, lon12;
	double ssig2, csig2, sbet2, cbet2, somg2, comg2, salp2, calp2, dn2;
	int outmask;

	lat2 = lon2 = azi2 = s12 = m12 = M12 = M21 = S12 = 0;
	/* Avoid warning about uninitialized B12. */
	B12 = AB1 = 0;
	outmask = 0;
	outmask |= plat2 ? GLatitude : GNone;
	outmask |= plon2 ? GLongitude : GNone;
	outmask |= pazi2 ? GAzimuth : GNone;
	outmask |= ps12 ? GDistance : GNone;
	outmask |= pm12 ? GReducedLength : GNone;
	outmask |= pM12 || pM21 ? GGeodesicScale : GNone;
	outmask |= pS12 ? GArea : GNone;

	outmask &= l->caps & OUT_ALL;
	if(!(flags & GArcMode || (l->caps & (GDistanceIn & OUT_ALL))))
	/* Impossible distance calculation requested */
		return NaN();

	if(flags & GArcMode){
		/* Interpret s12_a12 as spherical arc length */
		sig12 = s12_a12*DEG;
		sincosd(s12_a12, &ssig12, &csig12);
	}else{
		double tau12, s, c;

		/* Interpret s12_a12 as distance */
		tau12 = s12_a12 / (l->b * (1 + l->A1m1));
		s = sin(tau12);
		c = cos(tau12);
		/* tau2 = tau1 + tau12 */
		B12 = - SinCosSeries(1,
			l->stau1 * c + l->ctau1 * s,
			l->ctau1 * c - l->stau1 * s,
			l->C1pa, nC1p);
		sig12 = tau12 - (B12 - l->B11);
		ssig12 = sin(sig12);
		csig12 = cos(sig12);
		if(fabs(l->f) > 0.01){
			/* Reverted distance series is inaccurate for |f| > 1/100, so correct
			* sig12 with 1 Newton iteration.  The following table shows the
			* approximate maximum error for a = WGS_a() and various f relative to
			* GeodesicExact.
			*     erri = the error in the inverse solution (nm)
			*     errd = the error in the direct solution (series only) (nm)
			*     errda = the error in the direct solution (series + 1 Newton) (nm)
			*
			*       f     erri  errd errda
			*     -1/5    12e6 1.2e9  69e6
			*     -1/10  123e3  12e6 765e3
			*     -1/20   1110 108e3  7155
			*     -1/50  18.63 200.9 27.12
			*     -1/100 18.63 23.78 23.37
			*     -1/150 18.63 21.05 20.26
			*      1/150 22.35 24.73 25.83
			*      1/100 22.35 25.03 25.31
			*      1/50  29.80 231.9 30.44
			*      1/20   5376 146e3  10e3
			*      1/10  829e3  22e6 1.5e6
			*      1/5   157e6 3.8e9 280e6 */
			double serr;

			ssig2 = l->ssig1*csig12 + l->csig1*ssig12;
			csig2 = l->csig1*csig12 - l->ssig1*ssig12;
			B12 = SinCosSeries(1, ssig2, csig2, l->C1a, nC1);
			serr = (1 + l->A1m1)*(sig12 + (B12 - l->B11)) - s12_a12/l->b;
			sig12 = sig12 - serr/sqrt(1 + l->k2 * ssig2*ssig2);
			ssig12 = sin(sig12);
			csig12 = cos(sig12);
			/* Update B12 below */
		}
	}

	/* sig2 = sig1 + sig12 */
	ssig2 = l->ssig1*csig12 + l->csig1*ssig12;
	csig2 = l->csig1*csig12 - l->ssig1*ssig12;
	dn2 = sqrt(1 + l->k2 * ssig2*ssig2);
	if(outmask & (GDistance | GReducedLength | GGeodesicScale)){
		if(flags & GArcMode || fabs(l->f) > 0.01)
			B12 = SinCosSeries(1, ssig2, csig2, l->C1a, nC1);
		AB1 = (1 + l->A1m1)*(B12 - l->B11);
	}
	/* sin(bet2) = cos(alp0) * sin(sig2) */
	sbet2 = l->calp0 * ssig2;
	/* Alt: cbet2 = hypot(csig2, salp0 * ssig2); */
	cbet2 = hypot(l->salp0, l->calp0 * csig2);
	if(cbet2 == 0) /* I.e., salp0 = 0, csig2 = 0.  Break the degeneracy in this case */
		cbet2 = csig2 = tiny;
	/* tan(alp0) = cos(sig2)*tan(alp2) */
	salp2 = l->salp0; calp2 = l->calp0 * csig2; /* No need to normalize */

	if(outmask & GDistance)
		s12 = (flags & GArcMode) ?
			l->b*((1 + l->A1m1)*sig12 + AB1) :
			s12_a12;

	if(outmask & GLongitude){
		double E;

		E = copysign(1, l->salp0); /* east or west going? */
		/* tan(omg2) = sin(alp0) * tan(sig2) */
		somg2 = l->salp0 * ssig2;
		comg2 = csig2; /* No need to normalize */
		/* omg12 = omg2 - omg1 */
		omg12 = (flags & GUnrollLon) ?
			E*(sig12
			- (atan2(    ssig2, csig2) - atan2(    l->ssig1, l->csig1))
			+ (atan2(E * somg2, comg2) - atan2(E * l->somg1, l->comg1)))
			: atan2(somg2 * l->comg1 - comg2 * l->somg1,
				comg2 * l->comg1 + somg2 * l->somg1);
		lam12 = omg12 + l->A3c *
			(sig12 + (SinCosSeries(1, ssig2, csig2, l->C3a, nC3-1)
				- l->B31));
		lon12 = lam12/DEG;
		lon2 = (flags & GUnrollLon) ? l->lon1 + lon12 :
			AngNormalize(AngNormalize(l->lon1) + AngNormalize(lon12));
	}

	if(outmask & GLatitude)
		lat2 = atan2d(sbet2, l->f1*cbet2);

	if(outmask & GAzimuth)
		azi2 = atan2d(salp2, calp2);

	if(outmask & (GReducedLength | GGeodesicScale)){
		double B22, AB2, J12;

		B22 = SinCosSeries(1, ssig2, csig2, l->C2a, nC2);
		AB2 = (1 + l->A2m1)*(B22 - l->B21);
		J12 = (l->A1m1 - l->A2m1)*sig12 + (AB1 - AB2);
		if(outmask & GReducedLength)
			/* Add parens around (csig1 * ssig2) and (ssig1 * csig2) to ensure
			* accurate cancellation in the case of coincident points. */
			m12 = l->b * ((dn2 * (l->csig1 * ssig2) - l->dn1 * (l->ssig1 * csig2))
				- l->csig1 * csig2 * J12);
		if(outmask & GGeodesicScale){
			double t;

			t = l->k2*(ssig2 - l->ssig1)*(ssig2 + l->ssig1)/(l->dn1 + dn2);
			M12 = csig12 + (t*ssig2 -  csig2*J12) * l->ssig1/l->dn1;
			M21 = csig12 - (t*l->ssig1 - l->csig1*J12) * ssig2/dn2;
		}
	}

	if(outmask & GArea){
		double B42;
		double salp12, calp12;

		B42 = SinCosSeries(0, ssig2, csig2, l->C4a, nC4);
		if(l->calp0 == 0 || l->salp0 == 0){
			/* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
			salp12 = salp2*l->calp1 - calp2*l->salp1;
			calp12 = calp2*l->calp1 + salp2*l->salp1;
		}else{
			/* 
			 * tan(alp) = tan(alp0) * sec(sig)
			 * tan(alp2-alp1) = (tan(alp2) -tan(alp1)) / (tan(alp2)*tan(alp1)+1)
			 * = calp0 * salp0 * (csig1-csig2) / (salp0^2 + calp0^2 * csig1*csig2)
			 * If csig12 > 0, write
			 *   csig1 - csig2 = ssig12 * (csig1 * ssig12 / (1 + csig12) + ssig1)
			 * else
			 *   csig1 - csig2 = csig1 * (1 - csig12) + ssig12 * ssig1
			 * No need to normalize
			 */
			salp12 = l->calp0*l->salp0 *
				(csig12 <= 0 ? l->csig1*(1 - csig12) + ssig12*l->ssig1 :
				ssig12*(l->csig1*ssig12/(1 + csig12) + l->ssig1));
			calp12 = l->salp0*l->salp0 + l->calp0*l->calp0 * l->csig1*csig2;
		}
		S12 = l->c2*atan2(salp12, calp12) + l->A4*(B42 - l->B41);
	}

	if((outmask & GLatitude) && plat2)
		*plat2 = lat2;
	if((outmask & GLongitude) && plon2)
		*plon2 = lon2;
	if((outmask & GAzimuth) && pazi2)
		*pazi2 = azi2;
	if((outmask & GDistance) && ps12)
		*ps12 = s12;
	if((outmask & GReducedLength) && pm12)
		*pm12 = m12;
	if(outmask & GGeodesicScale){
		if(pM12)
			*pM12 = M12;
		if(pM21)
			*pM21 = M21;
	}
	if((outmask & GArea) && pS12)
		*pS12 = S12;

	return (flags & GArcMode) ? s12_a12 : sig12/DEG;
}

double
gendirectgeod(Geodesic *g,
	double lat1, double lon1, double azi1,
	uint flags, double s12_a12,
	double *lat2, double *lon2, double *azi2,
	double *s12, double *m12, double *M12, double *M21,
	double *S12)
{
	Geodline l;
	uint outmask;

	outmask = 0;
	outmask |= lat2 ? GLatitude : GNone;
	outmask |= lon2 ? GLongitude : GNone;
	outmask |= azi2 ? GAzimuth : GNone;
	outmask |= s12 ? GDistance : GNone;
	outmask |= m12 ? GReducedLength : GNone;
	outmask |= M12 ? GGeodesicScale : GNone;
	outmask |= S12 ? GArea : GNone;
	outmask |= flags&GArcMode ? GNone : GDistanceIn;
	initgeodline(&l, g, lat1, lon1, azi1, outmask);
	return geod_genposition(&l, flags, s12_a12, lat2, lon2, azi2, s12, m12, M12, M21, S12);
}

void
directgeod(Geodesic *g,
	double lat1, double lon1, double azi1,
	double s12,
	double *lat2, double *lon2, double *azi2)
{
	gendirectgeod(g, lat1, lon1, azi1, GNoflags, s12, lat2, lon2, azi2,
			nil, nil, nil, nil, nil);
}

static double
geninversegeod_int(Geodesic *g,
	double lat1, double lon1, double lat2, double lon2,
	double *ps12,
	double *psalp1, double *pcalp1,
	double *psalp2, double *pcalp2,
	double *pm12, double *pM12, double *pM21,
	double *pS12)
{
	double s12, m12, M12, M21, S12;
	double lon12, lon12s;
	double sbet1, cbet1, sbet2, cbet2, s12x, m12x;
	double dn1, dn2, lam12, slam12, clam12;
	double a12, sig12, calp1, salp1, calp2, salp2;
	double Ca[nC];
	double omg12, somg12, comg12;
	int latsign, lonsign, swapp;
	int meridian;
	uint outmask;

	s12 = m12 = M12 = M21 = S12 = 0;
	s12x = m12x = 0;
	a12 = sig12 = calp1 = salp1 = calp2 = salp2 = 0;
	/* somg12 > 1 marks that it needs to be calculated */
	omg12 = comg12 = 0;
	somg12 = 2;

	outmask = 0;
	outmask |= ps12 ? GDistance : GNone;
	outmask |= pm12 ? GReducedLength : GNone;
	outmask |= pM12 || pM21 ? GGeodesicScale : GNone;
	outmask |= pS12 ? GArea : GNone;

	outmask &= OUT_ALL;
	/*
	 * Compute longitude difference (AngDiff does this carefully).  Result is
	 * in [-180, 180] but -180 is only for west-going geodesics.  180 is for
	 * east-going and meridional geodesics.
	 */
	lon12 = AngDiff(lon1, lon2, &lon12s);
	/* Make longitude difference positive. */
	lonsign = lon12 >= 0 ? 1 : -1;
	/* If very close to being on the same half-meridian, then make it so. */
	lon12 = lonsign * AngRound(lon12);
	lon12s = AngRound((180 - lon12) - lonsign * lon12s);
	lam12 = lon12*DEG;
	if(lon12 > 90){
		sincosd(lon12s, &slam12, &clam12);
		clam12 = -clam12;
	}else
		sincosd(lon12, &slam12, &clam12);

	/* If really close to the equator, treat as on equator. */
	lat1 = AngRound(LatFix(lat1));
	lat2 = AngRound(LatFix(lat2));
	/*
	 * Swap points so that point with higher (abs) latitude is point 1
	 * If one latitude is a nan, then it becomes lat1.
	 */
	swapp = fabs(lat1) < fabs(lat2) ? -1 : 1;
	if(swapp < 0){
		lonsign *= -1;
		swap(&lat1, &lat2);
	}
	/* Make lat1 <= 0 */
	latsign = lat1 < 0 ? 1 : -1;
	lat1 *= latsign;
	lat2 *= latsign;
	/* Now we have
	*
	*     0 <= lon12 <= 180
	*     -90 <= lat1 <= 0
	*     lat1 <= lat2 <= -lat1
	*
	* longsign, swapp, latsign register the transformation to bring the
	* coordinates to this canonical form.  In all cases, 1 means no change was
	* made.  We make these transformations so that there are few cases to
	* check, e.g., on verifying quadrants in atan2.  In addition, this
	* enforces some symmetries in the results returned. */

	sincosd(lat1, &sbet1, &cbet1);
	sbet1 *= g->f1;
	/* Ensure cbet1 = +epsilon at poles */
	norm2(&sbet1, &cbet1);
	cbet1 = max(tiny, cbet1);

	sincosd(lat2, &sbet2, &cbet2);
	sbet2 *= g->f1;
	/* Ensure cbet2 = +epsilon at poles */
	norm2(&sbet2, &cbet2);
	cbet2 = max(tiny, cbet2);

	/*
	 * If cbet1 < -sbet1, then cbet2 - cbet1 is a sensitive measure of the
	 * |bet1| - |bet2|.  Alternatively (cbet1 >= -sbet1), abs(sbet2) + sbet1 is
	 * a better measure.  This logic is used in assigning calp2 in Lambda12.
	 * Sometimes these quantities vanish and in that case we force bet2 = +/-
	 * bet1 exactly.  An example where is is necessary is the inverse problem
	 * 48.522876735459 0 -48.52287673545898293 179.599720456223079643
	 * which failed with Visual Studio 10 (Release and Debug)
	 */
	if(cbet1 < -sbet1){
		if(cbet2 == cbet1)
			sbet2 = sbet2 < 0 ? sbet1 : -sbet1;
	}else{
		if (fabs(sbet2) == -sbet1)
			cbet2 = cbet1;
	}

	dn1 = sqrt(1 + g->ep2 * sbet1*sbet1);
	dn2 = sqrt(1 + g->ep2 * sbet2*sbet2);

	meridian = lat1 == -90 || slam12 == 0;

	if(meridian){
		/*
		 * Endpoints are on a single full meridian, so the geodesic might lie on
		 * a meridian.
		 */
	
		double ssig1, csig1, ssig2, csig2;
		/* Head to the target longitude */
		calp1 = clam12;
		salp1 = slam12;
		/* At the target we're heading north */
		calp2 = 1;
		salp2 = 0;
	
		/* tan(bet) = tan(sig) * cos(alp) */
		ssig1 = sbet1;
		csig1 = calp1 * cbet1;
		ssig2 = sbet2;
		csig2 = calp2 * cbet2;
	
		/* sig12 = sig2 - sig1 */
		sig12 = atan2(max(0, csig1 * ssig2 - ssig1 * csig2), csig1*csig2 + ssig1*ssig2);
		Lengths(g, g->n, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
			cbet1, cbet2, &s12x, &m12x, nil,
			(outmask & GGeodesicScale) ? &M12 : nil,
			(outmask & GGeodesicScale) ? &M21 : nil, Ca);
		/*
		 * Add the check for sig12 since zero length geodesics might yield m12 <
		 * 0.  Test case was
		 *
		 *    echo 20.001 0 20.001 0 | GeodSolve -i
		 *
		 * In fact, we will have sig12 > pi/2 for meridional geodesic which is
		 * not a shortest path.
		 */
		if(sig12 < 1 || m12x >= 0){
			/* Need at least 2, to handle 90 0 90 180 */
			if(sig12 < 3*tiny)
				sig12 = m12x = s12x = 0;
			m12x *= g->b;
			s12x *= g->b;
			a12 = sig12/DEG;
		}else
			/* m12 < 0, i.e., prolate and too close to anti-podal */
			meridian = 0;
	}

	if(!meridian && sbet1 == 0 && (g->f <= 0 || lon12s >= g->f*180)){
		/* Geodesic runs along equator */
		calp1 = calp2 = 0; salp1 = salp2 = 1;
		s12x = g->a * lam12;
		sig12 = omg12 = lam12 / g->f1;
		m12x = g->b*sin(sig12);
		if (outmask & GGeodesicScale)
			M12 = M21 = cos(sig12);
		a12 = lon12 / g->f1;
	}else if(!meridian){
		/*
		 * Now point1 and point2 belong within a hemisphere bounded by a
		 * meridian and geodesic is neither meridional or equatorial.
		 */
		/* Figure a starting point for Newton's method */
		double dnm;

		dnm = 0;
		sig12 = InverseStart(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2,
			lam12, slam12, clam12,
			&salp1, &calp1, &salp2, &calp2, &dnm, Ca);

		if(sig12 >= 0){
			/* Short lines (InverseStart sets salp2, calp2, dnm) */
			s12x = sig12 * g->b * dnm;
			m12x = dnm*dnm * g->b * sin(sig12/dnm);
			if(outmask & GGeodesicScale)
				M12 = M21 = cos(sig12/dnm);
			a12 = sig12/DEG;
			omg12 = lam12/(g->f1*dnm);
		}else{
			/*
			 * Newton's method.  This is a straightforward solution of f(alp1) =
			 * lambda12(alp1) - lam12 = 0 with one wrinkle.  f(alp) has exactly one
			 * root in the interval (0, pi) and its derivative is positive at the
			 * root.  Thus f(alp) is positive for alp > alp1 and negative for alp <
			 * alp1.  During the course of the iteration, a range (alp1a, alp1b) is
			 * maintained which brackets the root and with each evaluation of
			 * f(alp) the range is shrunk, if possible.  Newton's method is
			 * restarted whenever the derivative of f is negative (because the new
			 * value of alp1 is then further from the solution) or if the new
			 * estimate of alp1 lies outside (0,π); in this case, the new starting
			 * guess is taken to be (alp1a + alp1b)/2.
			 */
			double ssig1, csig1, ssig2, csig2, eps, domg12;
			double salp1a, calp1a, salp1b, calp1b;
			uint numit;
			int tripn, tripb;

			ssig1 = csig1 = ssig2 = csig2 = eps = domg12 = 0;
			numit = 0;
			/* Bracketing range */
			salp1a = tiny;
			calp1a = 1;
			salp1b = tiny;
			calp1b = -1;
			tripn = 0;
			tripb = 0;
			for(; numit < maxit2; ++numit){			
				/*
				 * the WGS84 test set: mean = 1.47, sd = 1.25, max = 16
				 * WGS84 and random input: mean = 2.85, sd = 0.60
				 */
				double dv, v;

				dv = 0;
				v = Lambda12(g, sbet1, cbet1, dn1, sbet2, cbet2, dn2, salp1, calp1,
					slam12, clam12,
					&salp2, &calp2, &sig12, &ssig1, &csig1, &ssig2, &csig2,
					&eps, &domg12, numit < maxit1, &dv, Ca);
				/* 2 * tol0 is approximately 1 ulp for a number in [0, pi]. */
				/* Reversed test to allow escape with NaNs */
				if(tripb || !(fabs(v) >= (tripn ? 8 : 1) * tol0))
					break;
				/* Update bracketing values */
				if(v > 0 && (numit > maxit1 || calp1/salp1 > calp1b/salp1b)){
					salp1b = salp1;
					calp1b = calp1;
				}else if(v < 0 && (numit > maxit1 || calp1/salp1 < calp1a/salp1a)){
					salp1a = salp1;
					calp1a = calp1;
				}
				if(numit < maxit1 && dv > 0){
					double dalp1, sdalp1, cdalp1, nsalp1;
			
					dalp1 = -v/dv;
					sdalp1 = sin(dalp1);
					cdalp1 = cos(dalp1);
					nsalp1 = salp1*cdalp1 + calp1*sdalp1;
					if(nsalp1 > 0 && fabs(dalp1) < PI){
						calp1 = calp1 * cdalp1 - salp1 * sdalp1;
						salp1 = nsalp1;
						norm2(&salp1, &calp1);
						/*
						 * In some regimes we don't get quadratic convergence because
						 * slope -> 0.  So use convergence conditions based on epsilon
						 * instead of sqrt(epsilon).
						 */
						tripn = fabs(v) <= 16 * tol0;
						continue;
					}
				}
				/*
				 * Either dv was not positive or updated value was outside legal
				 * range.  Use the midpoint of the bracket as the next estimate.
				 * This mechanism is not needed for the WGS84 ellipsoid, but it does
				 * catch problems with more eccentric ellipsoids.  Its efficacy is
				 * such for the WGS84 test set with the starting guess set to alp1 =
				 * 90deg:
				 * the WGS84 test set: mean = 5.21, sd = 3.93, max = 24
				 * WGS84 and random input: mean = 4.74, sd = 0.99
				 */
				salp1 = (salp1a + salp1b)/2;
				calp1 = (calp1a + calp1b)/2;
				norm2(&salp1, &calp1);
				tripn = 0;
				tripb = (fabs(salp1a - salp1) + (calp1a - calp1) < tolb ||
					fabs(salp1 - salp1b) + (calp1 - calp1b) < tolb);
				
			}
			Lengths(g, eps, sig12, ssig1, csig1, dn1, ssig2, csig2, dn2,
				cbet1, cbet2, &s12x, &m12x, nil,
				(outmask & GGeodesicScale) ? &M12 : nil,
				(outmask & GGeodesicScale) ? &M21 : nil, Ca);
			m12x *= g->b;
			s12x *= g->b;
			a12 = sig12/DEG;
			if(outmask & GArea){
				double sdomg12, cdomg12;

				/* omg12 = lam12 - domg12 */
				sdomg12 = sin(domg12);
				cdomg12 = cos(domg12);
				somg12 = slam12*cdomg12 - clam12*sdomg12;
				comg12 = clam12*cdomg12 + slam12*sdomg12;
			}
		}
	}

	if(outmask & GDistance)
		s12 = 0 + s12x;	/* Convert -0 to 0 */

	if(outmask & GReducedLength)
		m12 = 0 + m12x;	/* Convert -0 to 0 */

	if(outmask & GArea){
		double salp0, calp0, alp12;

		/* From Lambda12: sin(alp1) * cos(bet1) = sin(alp0) */
		salp0 = salp1 * cbet1;
		calp0 = hypot(calp1, salp1 * sbet1);	/* calp0 > 0 */
		if(calp0 != 0 && salp0 != 0){
			double ssig1, csig1, ssig2, csig2, k2, eps, A4;
			double B41, B42;

			/* From Lambda12: tan(bet) = tan(sig) * cos(alp) */
			ssig1 = sbet1;
			csig1 = calp1 * cbet1;
			ssig2 = sbet2;
			csig2 = calp2 * cbet2;
			k2 = calp0*calp0 * g->ep2;
			eps = k2/(2*(1 + sqrt(1 + k2)) + k2);
			/* Multiplier = a² * e² * cos(alpha0) * sin(alpha0). */
			A4 = g->a*g->a * calp0 * salp0 * g->e2;
			norm2(&ssig1, &csig1);
			norm2(&ssig2, &csig2);
			C4f(g, eps, Ca);
			B41 = SinCosSeries(0, ssig1, csig1, Ca, nC4);
			B42 = SinCosSeries(0, ssig2, csig2, Ca, nC4);
			S12 = A4*(B42 - B41);
		}else
			/* Avoid problems with indeterminate sig1, sig2 on equator */
			S12 = 0;

		if(!meridian && somg12 > 1){
			somg12 = sin(omg12);
			comg12 = cos(omg12);
		}
	
		if(!meridian && comg12 > -0.7071 && /* Lon difference not too big */
			sbet2 - sbet1 < 1.75){ /* Lat difference not too big */
			/*
			 * Use tan(Gamma/2) = tan(omg12/2) *
			 * 	(tan(bet1/2)+tan(bet2/2))/(1+tan(bet1/2)*tan(bet2/2))
			 * with tan(x/2) = sin(x)/(1+cos(x))
			 */
			double domg12, dbet1, dbet2;

			domg12 = 1 + comg12;
			dbet1 = 1 + cbet1;
			dbet2 = 1 + cbet2;
			alp12 = 2*atan2(somg12 * (sbet1 * dbet2 + sbet2 * dbet1),
				domg12*(sbet1 * sbet2 + dbet1 * dbet2));
		}else{
			/* alp12 = alp2 - alp1, used in atan2 so no need to normalize */
			double salp12, calp12;

			salp12 = salp2*calp1 - calp2*salp1;
			calp12 = calp2*calp1 + salp2*salp1;
			/*
			 * The right thing appears to happen if alp1 = +/-180 and alp2 = 0, viz
			 * salp12 = -0 and alp12 = -180.  However this depends on the sign
			 * being attached to 0 correctly.  The following ensures the correct
			 * behavior.
			 */
			if(salp12 == 0 && calp12 < 0){
				salp12 = tiny * calp1;
				calp12 = -1;
			}
			alp12 = atan2(salp12, calp12);
		}
		S12 += g->c2 * alp12;
		S12 *= swapp * lonsign*latsign;
		/* Convert -0 to 0 */
		S12 += 0;
	}

	/* Convert calp, salp to azimuth accounting for lonsign, swapp, latsign. */
	if(swapp < 0){
		swap(&salp1, &salp2);
		swap(&calp1, &calp2);
		if(outmask & GGeodesicScale)
			swap(&M12, &M21);
	}

	salp1 *= swapp*lonsign;
	calp1 *= swapp*latsign;
	salp2 *= swapp*lonsign;
	calp2 *= swapp*latsign;

	if(psalp1)
		*psalp1 = salp1;
	if(pcalp1)
		*pcalp1 = calp1;
	if(psalp2)
		*psalp2 = salp2;
	if(pcalp2)
		*pcalp2 = calp2;

	if (outmask & GDistance)
		*ps12 = s12;
	if (outmask & GReducedLength)
		*pm12 = m12;
	if(outmask & GGeodesicScale){
		if(pM12)
			*pM12 = M12;
		if(pM21)
			*pM21 = M21;
	}
	if(outmask & GArea)
		*pS12 = S12;
	/* Returned value in [0, 180] */
	return a12;
}


double
geninversegeod(Geodesic *g,
	double lat1, double lon1, double lat2, double lon2,
	double *ps12, double *pazi1, double *pazi2,
	double *pm12, double *pM12, double *pM21, double *pS12)
{
	double a12, salp1, calp1, salp2, calp2;

	a12 = geninversegeod_int(g, lat1, lon1, lat2, lon2, ps12,
		&salp1, &calp1, &salp2, &calp2,
		pm12, pM12, pM21, pS12);
	if(pazi1)
		*pazi1 = atan2d(salp1, calp1);
	if(pazi2)
		*pazi2 = atan2d(salp2, calp2);
	return a12;
}

void
inversegeod(Geodesic *g,
	double lat1, double lon1, double lat2, double lon2,
	double *s12,
	double *azi1, double *azi2)
{
	geninversegeod(g, lat1, lon1, lat2, lon2, s12, azi1, azi2, nil, nil, nil, nil);
}

void
initgeod(Geodesic *g, double a, double f)
{
	init();
	g->a = a;
	g->f = f;
	g->f1 = 1 - f;
	g->e2 = f*(2 - f);
	g->ep2 = g->e2 / g->f1*g->f1;   /* e2/(1 - e2) */
	g->n = f/(2 - f);
	g->b = a*g->f1;
	g->c2 = (a*a + g->b*g->b * (g->e2 == 0 ? 1 : (g->e2 > 0 ? atanh(sqrt(g->e2)) : atan(sqrt(-g->e2))) / sqrt(fabs(g->e2))))/2; /* authalic radius squared */
	/*
	 * The sig12 threshold for "really short".  Using the auxiliary sphere
	 * solution with dnm computed at (bet1 + bet2) / 2, the relative error in the
	 * azimuth consistency check is sig12^2 * abs(f) * min(1, 1-f/2) / 2.  (Error
	 * measured for 1/100 < b/a < 100 and abs(f) >= 1/1000.  For a given f and
	 * sig12, the max error occurs for lines near the pole.  If the old rule for
	 * computing dnm = (dn1 + dn2)/2 is used, then the error increases by a
	 * factor of 2.)  Setting this equal to epsilon gives sig12 = etol2.  Here
	 * 0.1 is a safety factor (error decreased by 100) and max(0.001, abs(f))
	 * stops etol2 getting too large in the nearly spherical case.
	 */
	g->etol2 = 0.1*tol2 / sqrt(max(0.001, fabs(f))*min(1, 1 - f/2)/2);
	A3coeff(g);
	C3coeff(g);
	C4coeff(g);
}
