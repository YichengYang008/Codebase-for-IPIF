struct Ran_FHDI{
//Description---------------------
// generate random double in [0, 1] 
// based on Numerical Recipes of Press et al. 2007
//--------------------------------
//example
//   Ran_FHDI myran(17);
//   myran.int64() returns a random 64-bit unsigned integer
//   myran.int32() returns a random 32-bit unsigned integer
//   myran.doub() returns a random double precision floating value in 0 and 1 
//--------------------------------
unsigned long long u, v, w; 
Ran_FHDI(unsigned long long j) : v(4101842887655102017LL), w(1) {
	u = j^v; int64();
	v = u; int64();
	w = v; int64(); 
}

//return 64-bit random integer
inline unsigned long long int64()
{
	u = u*2862933555777941757LL + 7046029254386353087LL;
	v ^= v >> 17;
	v ^= v << 31;
	v ^= v >> 8;
	w = 4294957665U*(w & 0xffffffff) + (w>>32);
	unsigned long long x = u^(u<<21);
	x^=x >> 35;
	x^=x << 4;
	return(x+v)^w;
}

//return random double-precision floating value in the range 0 to 1
inline double doub(){ return 5.42101086242752217E-20 * int64(); }

//return 32-bit random integer
inline unsigned int int32(){ return (unsigned int)int64(); }
};