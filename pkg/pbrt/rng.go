package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

const (
	PCG32DefaultState  = 0x853c49e6748fea9b
	PCG32DefaultStream = 0xda3e39cb94b95bdb
	PCG32Mult          = 0x5851f42d4c957f2d
)

type RandomNumberGenerator struct {
	state, inc uint64
}

func NewRandomNumberGenerator() *RandomNumberGenerator {
	return &RandomNumberGenerator{
		state: PCG32DefaultState,
		inc:   PCG32DefaultStream,
	}
}

func NewRNGWithSeed(seed uint64) *RandomNumberGenerator {
	r := NewRandomNumberGenerator()
	r.SetSequence(seed)
	return r
}

func (rng *RandomNumberGenerator) SetSequence(seed uint64) {
	rng.state = uint64(0)
	rng.inc = (seed << uint64(1)) | uint64(1)
	rng.UniformUInt32()
	rng.state += PCG32DefaultState
	rng.UniformUInt32()
}

func (rng *RandomNumberGenerator) UniformUInt32() uint32 {
	oldState := rng.state
	rng.state = oldState*PCG32Mult + rng.inc
	xorshifted := uint32((oldState>>uint(18) ^ oldState) >> uint(27))
	rot := uint32(oldState >> uint(59))
	return (xorshifted >> rot) | (xorshifted << ((rot + uint32(1)) & 31))
}

func (rng *RandomNumberGenerator) UniformUInt32B(b uint32) uint32 {
	threshold := (^b + 1) % b
	var r uint32
	for {
		r = rng.UniformUInt32()
		if r >= threshold {
			return r % b
		}
	}
}

func (rng *RandomNumberGenerator) UniformFloat() float64 {
	return math.Min(math.OneMinusEpsilon, float64(rng.UniformUInt32())*2.3283064365386963e-10)
}
