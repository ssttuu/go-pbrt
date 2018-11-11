package pbrt

type Sampler interface {
	GetSamplesPerPixel() int64
	Clone(seed uint64) Sampler
	Get1D() float64
	Get2D() *Point2f
	Get1DArray(n int) []float64
	Get2DArray(n int) []*Point2f
	GetCameraSample(pRaster *Point2i) *CameraSample
	StartPixel(p *Point2i)
	StartNextSample() bool
	RoundCount(n int) int
	Request2DArray(n int)
}

type sampler struct {
	SamplesPerPixel int64

	currentPixel                             *Point2i
	currentPixelSampleIndex                  int64
	samples1DArraySizes, samples2DArraySizes []int
	sampleArray1D                            [][]float64
	sampleArray2D                            [][]*Point2f

	array1DOffset, array2DOffset int
}

func NewSampler(samplesPerPixel int64) *sampler {
	sampleArray1D := make([][]float64, samplesPerPixel)
	for i := int64(0); i < samplesPerPixel; i++ {
		sampleArray1D[i] = make([]float64, samplesPerPixel)
	}

	sampleArray2D := make([][]*Point2f, samplesPerPixel)
	for i := int64(0); i < samplesPerPixel; i++ {
		sampleArray2D[i] = make([]*Point2f, samplesPerPixel)
	}

	return &sampler{
		SamplesPerPixel: samplesPerPixel,
		sampleArray1D:   sampleArray1D,
		sampleArray2D:   sampleArray2D,
	}
}

func (s *sampler) GetSamplesPerPixel() int64 {
	return s.SamplesPerPixel
}

func (s *sampler) Get1D() float64 {
	return 0
}

func (s *sampler) Get2D() *Point2f {
	return new(Point2f)
}

func (s *sampler) GetCameraSample(pRaster *Point2i) *CameraSample {
	return &CameraSample{
		pFilm: NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		time:  s.Get1D(),
		pLens: s.Get2D(),
	}
}

func (s *sampler) StartPixel(p *Point2i) {
	s.currentPixel = p
	s.currentPixelSampleIndex = 0
	// reset array offsets for next pixel sample
	s.array1DOffset = 0
	s.array2DOffset = 0
}

func (s *sampler) StartNextSample() bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	return s.currentPixelSampleIndex < s.SamplesPerPixel-1
}

func (s *sampler) SetSampleNumber(sampleNum int64) bool {
	s.array1DOffset = 0
	s.array2DOffset = 0
	s.currentPixelSampleIndex = sampleNum
	return s.currentPixelSampleIndex < s.SamplesPerPixel
}

func (s *sampler) Request1DArray(n int) {
	s.samples1DArraySizes = append(s.samples1DArraySizes, n)
	s.sampleArray1D = append(s.sampleArray1D, make([]float64, int64(n)*s.SamplesPerPixel, int64(n)*s.SamplesPerPixel))
}

func (s *sampler) Request2DArray(n int) {
	s.samples2DArraySizes = append(s.samples2DArraySizes, n)
	s.sampleArray2D = append(s.sampleArray2D, make([]*Point2f, int64(n)*s.SamplesPerPixel, int64(n)*s.SamplesPerPixel))
}

func (s *sampler) RoundCount(n int) int {
	return n
}

func (s *sampler) Get1DArray(n int) []float64 {
	//if s.array1DOffset == len(s.sampleArray1D) {
	//	return nil
	//}
	value := s.sampleArray1D[s.array1DOffset][s.currentPixelSampleIndex*int64(n):]
	s.array1DOffset++
	return value
}

func (s *sampler) Get2DArray(n int) []*Point2f {
	//if s.array2DOffset == len(s.sampleArray2D) {
	//	return nil
	//}
	value := s.sampleArray2D[s.array2DOffset][s.currentPixelSampleIndex*int64(n):]
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

func NewPixelSampler(samplesPerPixel int64, nSampledDimensions int) *PixelSampler {
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

func (s *PixelSampler) SetSampleNumber(sampleNum int64) bool {
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
		v = &Point2f{s.rng.UniformFloat(), s.rng.UniformFloat()}
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

func NewRandomSampler(ns, seed int) *RandomSampler {
	return &RandomSampler{
		sampler: NewSampler(int64(ns)),
		rng:     NewRNGWithSeed(uint64(seed)),
	}
}

func (s *RandomSampler) Get1D() float64 {
	return s.rng.UniformFloat()
}

func (s *RandomSampler) Get2D() *Point2f {
	return &Point2f{s.rng.UniformFloat(), s.rng.UniformFloat()}
}

func (s *RandomSampler) GetCameraSample(pRaster *Point2i) *CameraSample {
	return &CameraSample{
		pFilm: NewPoint2fFromPoint2i(pRaster).Add(s.Get2D()),
		time:  s.Get1D(),
		pLens: s.Get2D(),
	}
}

func (s *RandomSampler) Clone(seed uint64) Sampler {
	rs := &RandomSampler{
		sampler: (*sampler)(s.sampler),
		rng:     (*rng)(s.rng),
	}
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
			s.sampleArray2D[i][j] = &Point2f{s.rng.UniformFloat(), s.rng.UniformFloat()}
		}
	}

	s.sampler.StartPixel(p)
}
