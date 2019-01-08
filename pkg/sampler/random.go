package sampler

import "github.com/ssttuu/go-pbrt/pkg/pbrt"

// TODO: move out of core pbrt
type RandomSampler struct {
	*Sampler

	rng *pbrt.RandomNumberGenerator
}

func NewRandomSampler(ns int32, seed uint64) *RandomSampler {
	return &RandomSampler{
		Sampler: &Sampler{
			SamplesPerPixel: ns,
		},
		rng: pbrt.NewRNGWithSeed(seed),
	}
}

func (s *RandomSampler) Get1D() float64 {
	return s.rng.UniformFloat()
}

func (s *RandomSampler) Get2D() *pbrt.Point2f {
	return &pbrt.Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
}

func (s *RandomSampler) GetCameraSample(pRaster *pbrt.Point2i) *pbrt.CameraSample {
	return GetCameraSample(s, pRaster)
}

func (s *RandomSampler) Clone(seed uint64) pbrt.Sampler {
	rs := NewRandomSampler(s.SamplesPerPixel, seed)
	rs.rng.SetSequence(seed)
	return rs
}

func (s *RandomSampler) StartPixel(p *pbrt.Point2i) {
	for i := 0; i < len(s.sampleArray1D); i++ {
		for j := 0; j < len(s.sampleArray1D[i]); j++ {
			s.sampleArray1D[i][j] = s.rng.UniformFloat()
		}
	}

	for i := 0; i < len(s.sampleArray2D); i++ {
		for j := 0; j < len(s.sampleArray2D[i]); j++ {
			s.sampleArray2D[i][j] = pbrt.Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
		}
	}

	s.Sampler.StartPixel(p)
}
