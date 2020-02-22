#include <u.h>
#include <libc.h>
#include "geodesic.h"

typedef struct Geomodel Geomodel;

struct Geomodel
{
	double a;	/* equatorial radius */
	double f;	/* flattening */
};

Geomodel WGS84 = { 6378137, 1/298.257223563 };

void
usage(void)
{
	fprint(2, "usage: %s lat0 lon0 lat1 lon1\n", argv0);
	exits("usage");
}

void
main(int argc, char *argv[])
{
	Geodesic g;
	double lat0, lon0, lat1, lon1;
	double s12;

	ARGBEGIN{
	default: usage();
	}ARGEND;
	if(argc != 4)
		usage();
	lat0 = strtod(argv[0], nil);
	lon0 = strtod(argv[1], nil);
	lat1 = strtod(argv[2], nil);
	lon1 = strtod(argv[3], nil);
	initgeod(&g, WGS84.a, WGS84.f);
	inversegeod(&g, lat0, lon0, lat1, lon1, &s12, nil, nil);
	print("%gm\n", s12);
}
