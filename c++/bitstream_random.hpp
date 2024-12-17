//#######################################
//# Copyright (C) 2019-2020 Otmar Ertl. #
//# All rights reserved.                #
//#######################################

#ifndef _BIT_STREAM_RANDOM_HPP_
#define _BIT_STREAM_RANDOM_HPP_

#include "wyhash/wyhash.h"
#include "exponential_distribution.hpp"

#include <cmath>
#include <cassert>
#include <limits>
#include <memory>

static_assert(std::numeric_limits<double>::is_iec559, "Require std::numeric_limits<double>::is_iec559 to be true!");

static constexpr double maxInverse = 1. / (UINT64_C(1) << 53);
static constexpr double maxInverseHalf = 1. / (UINT64_C(1) << 54);

// uniform distributed double value from [0, 1) 
template<typename T> double getUniformDouble(T& bitstream) {
    return bitstream(53) * maxInverse;
}

// uniform distributed double value from [0, 0.5) 
template<typename T> double getUniformDoubleHalf(T& bitstream) {
    return bitstream(53) * maxInverseHalf;
}

template<typename T> double getExponential1(T& bitstream) {
    return -std::log1p(-getUniformDouble(bitstream));
}

template<typename T> double getGamma21(T& bitstream) {
    return getExponential1(bitstream) + getExponential1(bitstream);
}

template<typename T> bool getBernoulli(double successProbability, T& bitstream) {
    while(true) {
        if (successProbability == 0) return false;
        if (successProbability == 1) return true;
        bool b = successProbability > 0.5;
        if (bitstream()) return b;
        successProbability += successProbability;
        if (b) successProbability -= 1;
    }
}

template<typename T, typename V> bool getBernoulli(V successProbabilityNumerator, const V successProbabilityDenominator, T& bitstream) {
    assert(successProbabilityDenominator > 0);
    assert(successProbabilityNumerator <= successProbabilityDenominator);
    while(true) {
        if (successProbabilityNumerator == 0) return false;
        if (successProbabilityNumerator == successProbabilityDenominator) return true;
        bool b = successProbabilityNumerator > successProbabilityDenominator / 2;
        if (bitstream()) return b;
        successProbabilityNumerator += successProbabilityNumerator;
        if (b) successProbabilityNumerator -= successProbabilityDenominator;
    }
}

// returns smallest integer r such that v < 2^r
#if __GNUC__
static uint8_t getLog2Base(uint32_t v) {
    assert(v > 0);
    return static_cast<uint8_t>(std::numeric_limits<unsigned long>::digits - __builtin_clzl(v));
}
#else
static const uint8_t LogTable256[256] = {
    0, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8
};


static uint8_t getLog2Base(uint32_t v) {
    const uint32_t tt = v >> 16;
    if (tt != 0) {
        const uint32_t t = tt >> 8;
        if (t != 0) {
            assert(t <= 0xff);
            return 24 + LogTable256[t];
        }
        else {
            assert(tt <= 0xff);
            return 16 + LogTable256[tt];
        }
    } 
    else  {
        const uint32_t t = v >> 8;
        if (t != 0) {
            assert(t <= 0xff);
            return 8 + LogTable256[t];
        }
        else {
            assert(v <= 0xff);
            return LogTable256[v];
        }
    }
}
#endif

// A modified version of the method described in
// Lumbroso, Jeremie. "Optimal discrete uniform generation from coin flips, and applications." arXiv preprint arXiv:1304.1916 (2013).
// This method is very efficient in terms of consumed random bits.
template<typename T> uint32_t getUniformLumbroso(uint32_t n, T& bitstream) {
    assert(n >= 1);
    if (n == 1) {
        return 0;
    }

    uint8_t l = getLog2Base(n - 1); // l is in {1, 2,...,32}
    assert(l >= 1);
    uint64_t c = bitstream(l);

    if (c < n) {
        return static_cast<uint32_t>(c);
    }    
    uint64_t v = (UINT64_C(1) << l) - n;
    c -= n;

    while(true) {
        v <<= 1;
        c <<= 1;
        c += bitstream();
        if (v >= n) {
            if (c < n) {
                return static_cast<uint32_t>(c);
            }    
            v -= n;
            c -= n;
        }
    }
}

// Lemire, Daniel. "Fast random integer generation in an interval." ACM Transactions on Modeling and Computer Simulation (TOMACS) 29.1 (2019): 3.
// This method is very fast, but consumes more random bits than necessary.
template<typename T> uint32_t getUniformLemire(uint32_t s, T& bitstream) {
    uint64_t x = bitstream(32);
    uint64_t m = x * static_cast<uint64_t>(s);
    uint32_t l = static_cast<uint32_t>(m);
    if (l < s) {
        uint32_t t = -s % s;
        while (l < t) {
            x = bitstream(32);
            m = x * static_cast<uint64_t>(s);
            l = static_cast<uint32_t>(m);
        }
    }
    return static_cast<uint32_t>(m >> 32);
}

template<typename T> uint64_t getUniformPow2(uint8_t numBits, T& bitstream) {
    return bitstream(numBits);
}

class WyrandBitStream {

    uint64_t state;
    uint64_t hashBits;
    int availableBits;

    void nextHash() {
        availableBits += 64;
        hashBits = wyrand(&state);
    }
public:
    
    WyrandBitStream(const WyrandBitStream& p) = delete;
    WyrandBitStream& operator=(const WyrandBitStream&) = delete;
    WyrandBitStream(WyrandBitStream&& p) = default;
    WyrandBitStream& operator=(WyrandBitStream&&) = default;

    WyrandBitStream(uint64_t value, uint64_t seed) : state(wyhash64(value, seed)), hashBits(0), availableBits(0) {}
    
    WyrandBitStream(uint64_t value1, uint64_t value2, uint64_t seed) : hashBits(0), availableBits(0) {
        uint64_t data[2];
        data[0] = value1;
        data[1] = value2;
        state = wyhash(data, 2*sizeof(uint64_t), seed);
    }

    bool operator()() {
        assert(availableBits <= 63);
        if (availableBits == 0) {
            nextHash();
        }
        availableBits -= 1;
        bool result = hashBits >> availableBits;
        hashBits  &= ~(UINT64_C(0xFFFFFFFFFFFFFFFF) << availableBits);
        return result;
    }

    uint64_t operator()(uint8_t numBits) {
        assert(numBits >= 1);
        assert(numBits <= 64);
        assert(availableBits <= 63);
        uint64_t result = 0;
        availableBits -= numBits;
        if(availableBits < 0) {
            result |= hashBits << 1 << (- availableBits - 1);
            nextHash();
        }
        result |= hashBits >> availableBits;
        hashBits &= ~(UINT64_C(0xFFFFFFFFFFFFFFFF) << availableBits);
        return result;
    }
};

class TruncatedExponentialDistribution {
    double rate;
    double c1; // (exp(r) - 1) / r = expm1(r) / r
    double c2; // (-log(0.5*(1+exp(-r))) / r = -log1p(expm1(-r)*0.5) / r
    double c3; // expm1(r) / (exp(r) * r) = -expm1(-r) / r

public:

    TruncatedExponentialDistribution(double rate) : 
            rate(rate), 
            c1( (rate != 0) ? expm1(rate) / rate : 1),
            c2( (rate != 0) ? -log1p(expm1(-rate)*0.5) / rate : 0.5), 
            c3( (rate != 0) ? -expm1(-rate) / rate : 1)
    {
        assert(rate >= 0);
    }

    TruncatedExponentialDistribution() : TruncatedExponentialDistribution(0) {}

    template<typename T> 
    double operator()(T& bitstream) const {
        double x = getUniformDouble(bitstream) * c1;
        if (x < 1) {
            return x;
        }
        else {
            while(true) {
                x = getUniformDouble(bitstream);
                if (x <= c2) return x;
                double y = getUniformDoubleHalf(bitstream);
                if (y > 1 - x) { // 25% chance that this condition is satisfied 
                    x = 1 - x; 
                    y = 1 - y;
                }
                
                if ( x <= c3 * (1 - y)) return x;
                double c1y = c1 * y;
                if (c1y <= 1 - x) return x;
                if (c1y * rate <= expm1(rate * (1 - x))) return x;
            }
        }
    }
};

// based on Fisher-Yates shuffling
class PermutationStream {

    const uint32_t size;
	uint32_t idx;
	uint32_t versionCounter;

    const std::unique_ptr<std::pair<uint32_t, uint32_t>[]> permutationAndVersion;

public:

	PermutationStream(uint32_t _size) :
        size(_size), 
        idx(0), 
        versionCounter(0), 
        permutationAndVersion(new std::pair<uint32_t, uint32_t>[_size]) {}

	bool hasNext() const {
		return idx < size;
	}

	template <typename H>
    uint32_t next(H& hashBitStream) {
		const uint32_t k = idx + getUniformLemire(size - idx, hashBitStream);
        auto& permutationAndVersionK = permutationAndVersion[k];
        const auto& permutationAndVersionIdx = permutationAndVersion[idx];
        const uint32_t result = (permutationAndVersionK.second != versionCounter)?k:permutationAndVersionK.first;
        const uint32_t x = (permutationAndVersionIdx.second != versionCounter)?idx:permutationAndVersionIdx.first;
		permutationAndVersionK = std::make_pair(x, versionCounter);
		idx += 1;
        return result;
	}

    // must be called first before iterating over a new permutation using next-method
	void reset() {
		idx = 0;
		if (versionCounter == 0) {
            std::fill_n(permutationAndVersion.get(), size, std::make_pair(0,0));
		}
		versionCounter += 1;
	}
};

#endif // _BIT_STREAM_RANDOM_HPP_
