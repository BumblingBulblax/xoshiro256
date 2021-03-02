/*
 * xoshiro256.hpp
 *
 *  Created on: Feb 16, 2021
 *      Author: michael
 *
 *  This is based off of the C code provided on Sebastiano Vigna's website:
 *  http://prng.di.unimi.it/splitmix64.c
 *  http://prng.di.unimi.it/xoshiro256starstar.c
 *  http://prng.di.unimi.it/xoshiro256plus.c
 *
 *  This header file combines these scripts and makes them into classes.
 *  Additionally, there have been minor modifications to make this compatible
 *  with the probability distribution functions in <random>. However, I have found
 *  that these functions can be a bit slow, so I am adding some homemade functions
 *  for converting random variables of other distributions more efficiently. Xoshiro256**
 *  is used as the base class, as it is more precise, but Xoshiro256+ is explicitly
 *  generating floating point numbers, so it is also included via inheritance.
 *
 *  ----------------------Original SplitMix Comments----------------------
 *
 *  To the extent possible under law, the author has dedicated all copyright
 *  and related and neighboring rights to this software to the public domain
 *  worldwide. This software is distributed without any warranty.
 *
 *  See <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 *  This is a fixed-increment version of Java 8's SplittableRandom generator
 *  See http://dx.doi.org/10.1145/2714064.2660195 and
 *  http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html
 *
 *  It is a very fast generator passing BigCrush, and it can be useful if
 *  for some reason you absolutely want 64 bits of state.
 *
 *  ---------------------Original Xoshiro256** Comments---------------------
 *
 *  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 *  To the extent possible under law, the author has dedicated all copyright
 *  and related and neighboring rights to this software to the public domain
 *  worldwide. This software is distributed without any warranty.
 *
 *  See <http://creativecommons.org/publicdomain/zero/1.0/>
 *
 *  This is xoshiro256** 1.0, one of our all-purpose, rock-solid
 *  generators. It has excellent (sub-ns) speed, a state (256 bits) that is
 *  large enough for any parallel application, and it passes all tests we
 *  are aware of.
 *
 *  For generating just floating-point numbers, xoshiro256+ is even faster.
 *
 *  The state must be seeded so that it is not everywhere zero. If you have
 *  a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 *  output to fill s.
 *
 *  ---------------------Original Xoshiro256+ Comments---------------------
 *
 *  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
 *
 *	To the extent possible under law, the author has dedicated all copyright
 *	and related and neighboring rights to this software to the public domain
 *	worldwide. This software is distributed without any warranty.
 *
 *	See <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 *  This is xoshiro256+ 1.0, our best and fastest generator for floating-point
 *  numbers. We suggest to use its upper bits for floating-point
 *  generation, as it is slightly faster than xoshiro256++/xoshiro256**. It
 *  passes all tests we are aware of except for the lowest three bits,
 *  which might fail linearity tests (and just those), so if low linear
 *  complexity is not considered an issue (as it is usually the case) it
 *  can be used to generate 64-bit outputs, too.
 *
 *  We suggest to use a sign test to extract a random Boolean value, and
 *  right shifts to extract subsets of bits.
 *
 *  The state must be seeded so that it is not everywhere zero. If you have
 *  a 64-bit seed, we suggest to seed a splitmix64 generator and use its
 *  output to fill s.
 */
#include <cstdint>
#include <limits>
#include <chrono>
#include <cmath>
#include <sstream>
#include <iostream>
#ifndef XOSHIRO256_HPP_
#define XOSHIRO256_HPP_

/*
 * class declaration for 64-bit splitmix
 */
class splitmix64 {
public:
	uint64_t min() const; // returns 0
	uint64_t max() const; // returns the max uint64_t value
	splitmix64(uint64_t x0); // the constructor requires a seed
	uint64_t operator()(); // gets the next value. compatible with random's distributions
private:
	uint64_t x; // internal state
};

/*
 * class declaration for xoshiro256**
 */
class xoshiro256ss {
public:
	uint64_t min() const; // returns 0
	uint64_t max() const; // returns the max uint64_t value
	xoshiro256ss(); // default constructor with seeding from time and splitmix
	xoshiro256ss(uint64_t s0, uint64_t s1, uint64_t s2, uint64_t s3); // constructor with manual seeding
	virtual ~xoshiro256ss(){}; // destructor
	virtual uint64_t operator ()(); // gets the next value. compatible with random's distributions
	double uniform(double low, double high); // generates uniform reals in (a,b) using epsilon
	double exponential(double mean); // generates an exponential RV given the mean
	int geometric(double success); // generates a geometric RV... P(i failures) = p(1-p)^i
	void jump(); // this performs a jump
	void long_jump(); // this performs a larger jump
	uint64_t s[4]; // the state is four uint64_t
};

/*
 * class declaration for xoshiro256+... only the () operator is different
 */
class xoshiro256p : public xoshiro256ss {
public:
	uint64_t operator()() override; // gets the next value. compatible with random's distributions
	~xoshiro256p(){}; // destructor
};

/*
 * Rotation function using bit shifts. used in xoshiro256**
 */
inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

/*
 * splitmix64 constructor, requires a seed
 */
splitmix64::splitmix64(uint64_t x0){
	x=x0;
}

/*
 * get the next number from splitmix
 */
uint64_t splitmix64::operator ()() {
	uint64_t z = (x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

/*
 * splitmix min val
 */
uint64_t splitmix64::min() const{
	return 0;
}

/*
 * xoshiro max val
 */
uint64_t splitmix64::max() const{
	return std::numeric_limits<uint64_t>::max();
}

/*
 * xoshiro min val
 */
uint64_t xoshiro256ss::min() const{
	return 0;
}

/*
 * xoshiro max val
 */
uint64_t xoshiro256ss::max() const{
	return std::numeric_limits<uint64_t>::max();
}

/*
 * default xoshiro constructor. seeds a splitmix64 from the time, then
 * uses the first four outputs to seed xoshiro256**
 */
xoshiro256ss::xoshiro256ss(){
	splitmix64 seeder(std::chrono::high_resolution_clock::now()
									.time_since_epoch().count());
	s[0]=seeder();
	s[1]=seeder();
	s[2]=seeder();
	s[3]=seeder();
}

/*
 * specific xoshiro constructor, need to provide the four seeds
 */
xoshiro256ss::xoshiro256ss(uint64_t s0, uint64_t s1, uint64_t s2, uint64_t s3){
	s[0]=s0;
	s[1]=s1;
	s[2]=s2;
	s[3]=s3;
}

/*
 * get the next number from xoshiro256**
 */
uint64_t xoshiro256ss::operator()() {
	const uint64_t result = rotl(s[1] * 5, 7) * 9;

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);
	return result;
}

/*
 * get the next number from xoshiro256+
 */
uint64_t xoshiro256p::operator()() {
	const uint64_t result = s[0]+s[3];

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);
	return result;
}

/*
 * returns a uniform double in the open interval (low, high)
 */
double xoshiro256ss::uniform(double low, double high){
	// You could use epsilon to avoid n=0 or n=max, but it's faster to just check
	// and try again, if need be.
	uint64_t n = (*this)();
	while(n==0||n==std::numeric_limits<uint64_t>::max()){
		n = (*this)();
	}
	return low + (high-low)*n/((double)std::numeric_limits<uint64_t>::max());
}

/*
 * generates and exponential random variable with specified mean
 */
double xoshiro256ss::exponential(double mean){
	double r = uniform(0.0,1.0);
	return -mean*std::log(1-r);
}

/*
 * returns a geometric random variable (int)
 */
int xoshiro256ss::geometric(double success){
	double r = uniform(0.0,1.0);
	return std::ceil(-1+(std::log(1-r)/std::log(1-success)));
}

/*
 * jump function for xoshiro
 *
 * ----------------------Original Comments----------------------
 *
 * This is the jump function for the generator. It is equivalent
 * to 2^128 calls to next(); it can be used to generate 2^128
 * non-overlapping subsequences for parallel computations.
 */
void xoshiro256ss::jump() {
	const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			(*this)();
		}

	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

/*
 * long jump function for xoshiro
 *
 * ----------------------Original Comments----------------------
 *
 * This is the long-jump function for the generator. It is equivalent to
 * 2^192 calls to next(); it can be used to generate 2^64 starting points,
 * from each of which jump() will generate 2^64 non-overlapping
 * subsequences for parallel distributed computations.
 */

void xoshiro256ss::long_jump() {
	const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(unsigned int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			(*this)();
		}

	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

/*
 * converts uint64_t to strings. this is helpful for debugging.
 */
std::string UI64T2String(uint64_t input){
	uint64_t copy = input;
	std::string result = "";
	std::ostringstream ostrm;
	for(int i=0; i<64; i++){
		// Get the modulus wrt 2
		ostrm << (copy & 1);
		// append to the front of the result
		result = ostrm.str()+result;
		// shift down and clear the stream
		copy = copy >> 1;
		ostrm.str(std::string());
	}
	return result;
}

#endif /* XOSHIRO256_HPP_ */
