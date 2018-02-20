package geometry

type SpectrumGetter interface {
	NSamples() int
	GetSample(index int) float64
	GetSamples() []float64
}

type SpectrumSetter interface {
	SetSample(index int, value float64)
	SetSamples(values []float64)
}

type SpectrumAdder interface {
	Add(other SpectrumGetter)
}

type SpectrumSubtracter interface {
	Subtract(other SpectrumGetter)
}

type SpectrumMultiplier interface {
	Multiply(other SpectrumGetter)
}

type SpectrumDivider interface {
	Divide(other SpectrumGetter)
}

type SpectrumCloner interface {
	Clone() Spectrum
}

type Spectrum interface {
	SpectrumCloner
	SpectrumGetter
	SpectrumSetter
	SpectrumAdder
	SpectrumSubtracter
	SpectrumMultiplier
	SpectrumDivider
}

func SpectrumAdd(s1 SpectrumCloner, s2 SpectrumGetter) Spectrum {
	s3 := s1.Clone()
	s3.Add(s2)
	return s3
}

func SpectrumSubtract(s1 SpectrumCloner, s2 SpectrumGetter) Spectrum {
	s3 := s1.Clone()
	s3.Subtract(s2)
	return s3
}

func SpectrumMultiply(s1 SpectrumCloner, s2 SpectrumGetter) Spectrum {
	s3 := s1.Clone()
	s3.Multiply(s2)
	return s3
}

func SpectrumDivide(s1 SpectrumCloner, s2 SpectrumGetter) Spectrum {
	s3 := s1.Clone()
	s3.Divide(s2)
	return s3
}

func SpectrumEquals(s1, s2 SpectrumGetter) bool {
	for i := 0; i < s1.NSamples(); i++ {
		if s1.GetSample(i) != s2.GetSample(i) {
			return false
		}
	}
	return true
}

func SpectrumNotEquals(s1, s2 SpectrumGetter) bool {
	for i := 0; i < s1.NSamples(); i++ {
		if s1.GetSample(i) != s2.GetSample(i) {
			return true
		}
	}
	return false
}

type RGBSpectrum struct {
	samples []float64
}

func NewRGBSpectrum(samples []float64) RGBSpectrum {
	return RGBSpectrum{samples:samples}
}

func (s *RGBSpectrum) Clone() Spectrum {
	return &RGBSpectrum{samples: s.samples}
}

func (s *RGBSpectrum) NSamples() int {
	return len(s.samples)
}

func (s *RGBSpectrum) GetSample(index int) float64 {
	return s.samples[index]
}

func (s *RGBSpectrum) GetSamples() []float64 {
	return s.samples
}

func (s *RGBSpectrum) SetSample(index int, value float64) {
	s.samples[index] = value
}

func (s *RGBSpectrum) SetSamples(values []float64) {
	s.samples = values
}

func (s *RGBSpectrum) Add(other SpectrumGetter) {
	for i := 0; i < s.NSamples(); i++ {
		s.SetSample(i, s.GetSample(i) + other.GetSample(i))
	}
}

func (s *RGBSpectrum) Subtract(other SpectrumGetter) {
	for i := 0; i < s.NSamples(); i++ {
		s.SetSample(i, s.GetSample(i) - other.GetSample(i))
	}
}

func (s *RGBSpectrum) Multiply(other SpectrumGetter) {
	for i := 0; i < s.NSamples(); i++ {
		s.SetSample(i, s.GetSample(i) * other.GetSample(i))
	}
}

func (s *RGBSpectrum) Divide(other SpectrumGetter) {
	for i := 0; i < s.NSamples(); i++ {
		s.SetSample(i, s.GetSample(i) / other.GetSample(i))
	}
}

// TODO
type SampledSpectrum struct {
	samples []float64
}

//func (s *SampledSpectrum) Clone() Spectrum {
//	return &SampledSpectrum{samples:s.samples}
//}
