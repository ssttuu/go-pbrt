package sampler

import "github.com/ssttuu/go-pbrt/pkg/pbrt"

type Sampler struct {
	SamplesPerPixel int32

	currentPixel                             pbrt.Point2i
	currentPixelSampleIndex                  int32
	samples1DArraySizes, samples2DArraySizes []int32
	sampleArray1D                            [][]float64
	sampleArray2D                            [][]pbrt.Point2f

	array1DOffset, array2DOffset int
}

func (s *Sampler) GetSamplesPerPixel() int32 {
	return s.SamplesPerPixel
}

func (s *Sampler) StartPixel(p *pbrt.Point2i) {
	s.currentPixel = *p
	s.currentPixelSampleIndex = 0
	// reset array offsets for next pixel sample
	s.array1DOffset = 0
	s.array2DOffset = 0
}

func (s *Sampler) StartNextSample() bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	s.currentPixelSampleIndex += 1
	return s.currentPixelSampleIndex < s.SamplesPerPixel
}

func (s *Sampler) SetSampleNumber(sampleNum int32) bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	s.currentPixelSampleIndex = sampleNum
	return s.currentPixelSampleIndex < s.SamplesPerPixel
}

func (s *Sampler) Request1DArray(n int32) {
	s.samples1DArraySizes = append(s.samples1DArraySizes, n)
	s.sampleArray1D = append(s.sampleArray1D, make([]float64, n*s.SamplesPerPixel, n*s.SamplesPerPixel))
}

func (s *Sampler) Request2DArray(n int32) {
	s.samples2DArraySizes = append(s.samples2DArraySizes, n)
	s.sampleArray2D = append(s.sampleArray2D, make([]pbrt.Point2f, n*s.SamplesPerPixel, n*s.SamplesPerPixel))
}

func (s *Sampler) RoundCount(n int32) int32 {
	return n
}

func (s *Sampler) Get1DArray(n int32) []float64 {
	if s.array1DOffset == len(s.sampleArray1D) {
		return nil
	}
	value := s.sampleArray1D[s.array1DOffset][s.currentPixelSampleIndex*n:]
	s.array1DOffset++
	return value
}

func (s *Sampler) Get2DArray(n int32) []pbrt.Point2f {
	if s.array2DOffset == len(s.sampleArray2D) {
		return nil
	}
	value := s.sampleArray2D[s.array2DOffset][s.currentPixelSampleIndex*n:]
	s.array2DOffset++
	return value
}

func GetCameraSample(s pbrt.Sampler, pRaster *pbrt.Point2i) *pbrt.CameraSample {
	return pbrt.NewCameraSample(
		pbrt.NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		s.Get2D(),
		s.Get1D(),
	)
}