package pbrt

import "github.com/stupschwartz/go-pbrt/pkg/math"

const (
	PCG32DefaultState  = 0x853c49e6748fea9b
	PCG32DefaultStream = 0xda3e39cb94b95bdb
	PCG32Mult          = 0x5851f42d4c957f2d
)

type rng struct {
	state, inc uint64
}

func NewRNG() *rng {
	return &rng{
		state: PCG32DefaultState,
		inc:   PCG32DefaultStream,
	}
}

func NewRNGWithSeed(seed uint64) *rng {
	r := NewRNG()

	r.SetSequence(seed)

	return r
}

func (r *rng) SetSequence(seed uint64) {
	r.state = uint64(0)
	r.inc = (seed << uint64(1)) | uint64(1)
	r.UniformUInt32()
	r.state += PCG32DefaultState
	r.UniformUInt32()
}

func (r *rng) UniformUInt32() uint32 {
	oldState := r.state
	r.state = oldState*PCG32Mult + r.inc
	xorshifted := uint32((oldState>>uint(18) ^ oldState) >> uint(27))
	rot := uint32(oldState >> uint(59))
	return (xorshifted >> rot) | (xorshifted << ((rot + uint32(1)) & 31))
}

func (r *rng) UniformFloat() float64 {
	return math.Min(math.OneMinusEpsilon, float64(r.UniformUInt32())*2.3283064365386963e-10)
}
