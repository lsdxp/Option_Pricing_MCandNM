#ifndef SOBOL_H_150124
#define SOBOL_H_150124

#include <vector>
#include <cassert>

typedef unsigned int uint32;
typedef unsigned long int uint64;

class Sobol{
    /* ---Sobol' sequence generator---
     * 	Dimension s is defined by the size of primitive polynomial's coefficient.
     * 	In this class,
     *		vector<char> polynomial(2);
     *		polynomial[0] = 0;
     *		polynomial[1] = 1;
     *	means x^2 + 0*x^1 + 1*x^0.
     *	(Originally the last coefficient must be 1)
     * 	*/
public:
    Sobol();
    Sobol(const uint32 bit_size, const std::vector<char>& polynomial);
    // Initialize
    void init(const uint32 bit_size, const std::vector<char>& polynomial);
    // Generate i-th number
    uint64 get(uint64 i);
private:
    uint32 shift;	// effective bit size
    uint32 s;	// polynomial's dimension
    std::vector<char> polynomial;	// primitive polynomial for generating v_array
    std::vector<uint64> v_array;	// LUT of v
    
    // Precompute LUT of v (v_array)
    void precompute();
};

#endif