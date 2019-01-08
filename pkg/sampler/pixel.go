package sampler

import "github.com/ssttuu/go-pbrt/pkg/pbrt"

type PixelSampler struct {
	*Sampler

	nSampledDimensions                     int
	samples1D                              [][]float64
	samples2D                              [][]pbrt.Point2f
	current1DDimension, current2DDimension int
	rng                                    *pbrt.RandomNumberGenerator
}

func NewPixelSampler(samplesPerPixel int32, nSampledDimensions int) *PixelSampler {
	ps := &PixelSampler{
		Sampler: &Sampler{
			SamplesPerPixel: samplesPerPixel,
		},
		nSampledDimensions: nSampledDimensions,
		rng:                pbrt.NewRandomNumberGenerator(),
	}
	for i := 0; i < nSampledDimensions; i++ {
		ps.samples1D = append(ps.samples1D, make([]float64, samplesPerPixel, samplesPerPixel))
		ps.samples2D = append(ps.samples2D, make([]pbrt.Point2f, samplesPerPixel, samplesPerPixel))
	}
	return ps
}

func (s *PixelSampler) Clone(seed uint64) pbrt.Sampler {
	return s.clone(seed)
}

func (s *PixelSampler) clone(seed uint64) *PixelSampler {
	ps := NewPixelSampler(s.Sampler.SamplesPerPixel, s.nSampledDimensions)
	for i := 0; i < s.nSampledDimensions; i++ {
		copy(ps.samples1D[i], s.samples1D[i])
		copy(ps.samples2D[i], s.samples2D[i])
	}
	ps.rng.SetSequence(seed)
	return ps
}

func (s *PixelSampler) GetCameraSample(pRaster *pbrt.Point2i) *pbrt.CameraSample {
	return GetCameraSample(s, pRaster)
}

func (s *PixelSampler) StartNextSample() bool {
	s.current1DDimension = 0
	s.current2DDimension = 0
	return s.Sampler.StartNextSample()
}

func (s *PixelSampler) SetSampleNumber(sampleNum int32) bool {
	s.current1DDimension = 0
	s.current2DDimension = 0
	return s.Sampler.SetSampleNumber(sampleNum)
}

func (s *PixelSampler) Get1D() float64 {
	var v float64
	if s.current1DDimension < len(s.samples1D) {
		v = s.samples1D[s.current1DDimension][s.currentPixelSampleIndex]
		s.current1DDimension++
	} else {
		v = s.rng.UniformFloat()
	}
	return v
}

func (s *PixelSampler) Get2D() *pbrt.Point2f {
	var v *pbrt.Point2f
	if s.current2DDimension < len(s.samples2D) {
		v = &s.samples2D[s.current2DDimension][s.currentPixelSampleIndex]
		s.current2DDimension++
	} else {
		v = &pbrt.Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
	}
	return v
}
