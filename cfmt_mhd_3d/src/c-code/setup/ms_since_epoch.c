#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>

/* Return number of ms since 1.1.1970, in a 64 bit integer.
 * (with 2^64 ms ~ 6 * 10^8 years, this should be sufficiently overflow safe)
 */
void ms_since_epoch(int64_t *ms) {
	struct timeval tv;
	const int64_t million = 1000000;
	if (gettimeofday(&tv, NULL) != 0) {
		perror("gettimeofday");
		exit(1);
	}
	*ms = (int64_t) (tv.tv_sec) * million + (int64_t)(tv.tv_usec);
}
