package sampler

import "github.com/ssttuu/go-pbrt/pkg/pbrt"

type Stratified struct {
	*PixelSampler

	xPixelSamples, yPixelSamples int32
	jitterSamples                bool
}

func NewStratified(xSamples, ySamples int32, jitterSamples bool, nSampledDimensions int) pbrt.Sampler {
	return &Stratified{
		PixelSampler: NewPixelSampler(xSamples * ySamples, nSampledDimensions),
		xPixelSamples: xSamples,
		yPixelSamples: ySamples,
		jitterSamples: jitterSamples,
	}
}

func (s *Stratified) StartPixel(p *pbrt.Point2i) {
	// Generate single stratified samples for the pixel
	for i := 0; i < len(s.PixelSampler.samples1D); i++ {
		pbrt.StratifiedSample1D(s.PixelSampler.samples1D[i], s.xPixelSamples * s.yPixelSamples, s.rng, s.jitterSamples)
		pbrt.ShuffleSamples1D(s.samples1D[i], s.xPixelSamples * s.yPixelSamples, 1, s.rng)
	}
	for i := 0; i < len(s.PixelSampler.samples2D); i++ {
		pbrt.StratifiedSample2D(s.PixelSampler.samples2D[i], s.xPixelSamples, s.yPixelSamples, s.rng, s.jitterSamples)
		pbrt.ShuffleSamples2D(s.samples2D[i], s.xPixelSamples * s.yPixelSamples, 1, s.rng)
	}

	// Generate arrays of stratified samples for the pixel
	for i := 0; i < len(s.samples1DArraySizes); i++ {
		for j := int32(0); j < s.SamplesPerPixel; j++ {
			count := s.samples1DArraySizes[i]
			pbrt.StratifiedSample1D(s.sampleArray1D[i][j * count:], count, s.rng, s.jitterSamples)
			pbrt.ShuffleSamples1D(s.sampleArray1D[i][j * count:], count, 1, s.rng)
		}
	}
	for i := 0; i < len(s.samples2DArraySizes); i++ {
		for j := int32(0); j < s.SamplesPerPixel; j++ {
			count := s.samples2DArraySizes[i]
			pbrt.LatinHypercube2D(s.sampleArray2D[i][j*count:count], s.rng)
		}
	}

	s.PixelSampler.StartPixel(p)
}

func (s *Stratified) Clone(seed uint64) pbrt.Sampler {
	return &Stratified{
		PixelSampler: s.PixelSampler.clone(seed),
		xPixelSamples: s.xPixelSamples,
		yPixelSamples: s.yPixelSamples,
		jitterSamples: s.jitterSamples,
	}
}