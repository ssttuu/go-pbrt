//go:generate mockgen -source=sampler.go -destination=sampler.mock.go -package=pbrt

package pbrt

type Sampler interface {
	GetSamplesPerPixel() int32
	Clone(seed uint64) Sampler
	Get1D() float64
	Get2D() *Point2f
	Get1DArray(n int32) []float64
	Get2DArray(n int32) []Point2f
	GetCameraSample(pRaster *Point2i) *CameraSample
	StartPixel(p *Point2i)
	StartNextSample() bool
	RoundCount(n int32) int32
	Request2DArray(n int32)
}

type sampler struct {
	SamplesPerPixel int32

	currentPixel                             Point2i
	currentPixelSampleIndex                  int32
	samples1DArraySizes, samples2DArraySizes []int32
	sampleArray1D                            [][]float64
	sampleArray2D                            [][]Point2f

	array1DOffset, array2DOffset int
}

func (s *sampler) GetSamplesPerPixel() int32 {
	return s.SamplesPerPixel
}

func GetCameraSample(s Sampler, pRaster *Point2i) *CameraSample {
	return &CameraSample{
		pFilm: NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		time:  s.Get1D(),
		pLens: s.Get2D(),
	}
}

func (s *sampler) StartPixel(p *Point2i) {
	s.currentPixel = *p
	s.currentPixelSampleIndex = 0
	// reset array offsets for next pixel sample
	s.array1DOffset = 0
	s.array2DOffset = 0
}

func (s *sampler) StartNextSample() bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	s.currentPixelSampleIndex += 1
	return s.currentPixelSampleIndex < s.SamplesPerPixel
}

func (s *sampler) SetSampleNumber(sampleNum int32) bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	s.currentPixelSampleIndex = sampleNum
	return s.currentPixelSampleIndex < s.SamplesPerPixel
}

func (s *sampler) Request1DArray(n int32) {
	s.samples1DArraySizes = append(s.samples1DArraySizes, n)
	s.sampleArray1D = append(s.sampleArray1D, make([]float64, n*s.SamplesPerPixel, n*s.SamplesPerPixel))
}

func (s *sampler) Request2DArray(n int32) {
	s.samples2DArraySizes = append(s.samples2DArraySizes, n)
	s.sampleArray2D = append(s.sampleArray2D, make([]Point2f, n*s.SamplesPerPixel, n*s.SamplesPerPixel))
}

func (s *sampler) RoundCount(n int32) int32 {
	return n
}

func (s *sampler) Get1DArray(n int32) []float64 {
	if s.array1DOffset == len(s.sampleArray1D) {
		return nil
	}
	value := s.sampleArray1D[s.array1DOffset][s.currentPixelSampleIndex*n:]
	s.array1DOffset++
	return value
}

func (s *sampler) Get2DArray(n int32) []Point2f {
	if s.array2DOffset == len(s.sampleArray2D) {
		return nil
	}
	value := s.sampleArray2D[s.array2DOffset][s.currentPixelSampleIndex*n:]
	s.array2DOffset++
	return value
}

type PixelSampler struct {
	*sampler

	samples1D                              [][]float64
	samples2D                              [][]*Point2f
	current1DDimension, current2DDimension int
	rng                                    *rng
}

func NewPixelSampler(samplesPerPixel int32, nSampledDimensions int) *PixelSampler {
	ps := &PixelSampler{
		sampler: &sampler{
			SamplesPerPixel: samplesPerPixel,
		},
		rng: NewRNG(),
	}
	for i := 0; i < nSampledDimensions; i++ {
		ps.samples1D = append(ps.samples1D, make([]float64, samplesPerPixel, samplesPerPixel))
		ps.samples2D = append(ps.samples2D, make([]*Point2f, samplesPerPixel, samplesPerPixel))
	}
	return ps
}

func (s *PixelSampler) StartNextSample() bool {
	s.current1DDimension = 0
	s.current2DDimension = 0
	return s.sampler.StartNextSample()
}

func (s *PixelSampler) SetSampleNumber(sampleNum int32) bool {
	s.current1DDimension = 0
	s.current2DDimension = 0
	return s.sampler.SetSampleNumber(sampleNum)
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

func (s *PixelSampler) Get2D() *Point2f {
	var v *Point2f
	if s.current2DDimension < len(s.samples2D) {
		v = s.samples2D[s.current2DDimension][s.currentPixelSampleIndex]
		s.current2DDimension++
	} else {
		v = &Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
	}
	return v
}

type GlobalSampler struct {
	*sampler
}

// TODO: move out of core pbrt
type RandomSampler struct {
	*sampler

	rng *rng
}

func NewRandomSampler(ns int32, seed uint64) *RandomSampler {
	return &RandomSampler{
		sampler: &sampler{
			SamplesPerPixel: ns,
		},
		rng: NewRNGWithSeed(seed),
	}
}

func (s *RandomSampler) Get1D() float64 {
	return s.rng.UniformFloat()
}

func (s *RandomSampler) Get2D() *Point2f {
	return &Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
}

func (s *RandomSampler) GetCameraSample(pRaster *Point2i) *CameraSample {
	return &CameraSample{
		pFilm: NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		time:  s.Get1D(),
		pLens: s.Get2D(),
	}
}

func (s *RandomSampler) Clone(seed uint64) Sampler {
	rs := NewRandomSampler(s.SamplesPerPixel, seed)
	rs.rng.SetSequence(seed)
	return rs
}

func (s *RandomSampler) StartPixel(p *Point2i) {
	for i := 0; i < len(s.sampleArray1D); i++ {
		for j := 0; j < len(s.sampleArray1D[i]); j++ {
			s.sampleArray1D[i][j] = s.rng.UniformFloat()
		}
	}

	for i := 0; i < len(s.sampleArray2D); i++ {
		for j := 0; j < len(s.sampleArray2D[i]); j++ {
			s.sampleArray2D[i][j] = Point2f{X: s.rng.UniformFloat(), Y: s.rng.UniformFloat()}
		}
	}

	s.sampler.StartPixel(p)
}
