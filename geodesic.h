/*
 * This an implementation in C of the geodesic algorithms described in
 * C. F. F. Karney, Algorithms for geodesics, J. Geodesy 87, pp. 43-55 (2013).
 *
 * The principal advantages of these algorithms over previous ones (e.g.,
 * Vincenty, 1975) are
 * - accurate to round off for |f| < 1/50
 * - the solution of the inverse problem is always found
 * - differential and integral properties of geodesics are computed
 *
 * The shortest path between two points on the ellipsoid at lat1, lon1
 * and lat2, lon2 is called the geodesic.  Its length is
 * s12 and the geodesic from point 1 to point 2 has forward azimuths
 * azi1 and azi2 at the two end points.
 *
 * Ported from PROJ v6.2.1 to Plan 9 on 22dec2019.
 */
#define DEG 0.0174532925199

enum {
	GeodesicOrder	= 6,
	nA1		= GeodesicOrder,
	nC1		= GeodesicOrder,
	nC1p		= GeodesicOrder,
	nA2		= GeodesicOrder,
	nC2		= GeodesicOrder,
	nA3		= GeodesicOrder,
	nA3x		= nA3,
	nC3		= GeodesicOrder,
	nC3x		= nC3*(nC3 - 1) / 2,
	nC4		= GeodesicOrder,
	nC4x		= nC4*(nC4 + 1) / 2,
	nC		= GeodesicOrder+1,
};

enum {
	GNone		= 0,
	GLatitude	= 1<<0,
	GLongitude	= 1<<1,
	GAzimuth	= 1<<2,
	GDistance	= 1<<3,
	GDistanceIn	= 1<<4,	/* Allow distance as input  */
	GReducedLength	= 1<<5,
	GGeodesicScale	= 1<<6,
	GArea		= 1<<7,
	GAll		= 0xFF,

	GNoflags	= 0,
	GArcMode	= 1<<0,
	GUnrollLon	= 1<<1
};

enum {
	CapNone	= 0,
	CapC1	= 1<<0,
	CapC1p	= 1<<1,
	CapC2	= 1<<2,
	CapC3	= 1<<3,
	CapC4	= 1<<4,
	CapAll	= 0x1F,
	OUT_ALL	= 0xFF
};

typedef struct Geodesic Geodesic;
typedef struct Geodline Geodline;

struct Geodesic
{
	double a;	/* the equatorial radius */
	double f;	/* the flattening */
	double f1;
	double e2, ep2;
	double n, b, c2, etol2;
	double A3x[6];
	double C3x[15];
	double C4x[21];
};

struct Geodline
{
	double lat1;	/* the starting latitude */
	double lon1;	/* the starting longitude */
	double azi1;	/* the starting azimuth */
	double a;	/* the equatorial radius */
	double f;	/* the flattening */
	double salp1;	/* sine of azi1 */
	double calp1;	/* cosine of azi1 */
	double a13;	/* arc length to reference point */
	double s13;	/* distance to reference point */
	double b, c2, f1, salp0, calp0, k2,
		ssig1, csig1, dn1, stau1, ctau1, somg1, comg1,
		A1m1, A2m1, A3c, B11, B21, B31, A4, B41;
	double C1a[6+1];
	double C1pa[6+1];
	double C2a[6+1];
	double C3a[6];
	double C4a[6];
	uint caps;	/* the capabilities */
};


void initgeod(Geodesic*, double, double);
void directgeod(Geodesic*, double, double, double, double, double*, double*, double*);
double gendirectgeod(Geodesic*, double, double, double, uint, double, double*,
	double*, double*, double*, double*, double*, double*, double*);
void inversegeod(Geodesic*, double, double, double, double, double*, double*, double*);
double geninversegeod(Geodesic*, double, double, double, double, double*,
	double*, double*, double*, double*, double*, double*);
void initgeodline(Geodline*, Geodesic*, double, double, double, uint);
double geod_genposition(Geodline*, uint, double, double*, double*, double*,
	double*, double*, double*, double*, double*);
