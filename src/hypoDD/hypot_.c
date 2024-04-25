#ifndef lint
#endif /* lint */

#include <math.h>

/* f77-callable interface to hypot function */
double
hypot_(a, b)
	float	*a, *b;
{
	return hypot((double)*a, (double)*b);
}
