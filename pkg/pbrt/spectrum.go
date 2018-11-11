package pbrt

import (
	"math"
	"fmt"
)

const (
	nSpectralSamples  = 60
	nRGB2SpectSamples = 32
)

//const (
//	X, Y, Z                                     *SampledSpectrum
//	rgbRefl2SpectWhite, rgbRefl2SpectCyan       *SampledSpectrum
//	rgbRefl2SpectMagenta, rgbRefl2SpectYellow   *SampledSpectrum
//	rgbRefl2SpectRed, rgbRefl2SpectGreen        *SampledSpectrum
//	rgbRefl2SpectBlue                           *SampledSpectrum
//	rgbIllum2SpectWhite, rgbIllum2SpectCyan     *SampledSpectrum
//	rgbIllum2SpectMagenta, rgbIllum2SpectYellow *SampledSpectrum
//	rgbIllum2SpectRed, rgbIllum2SpectGreen      *SampledSpectrum
//	rgbIllum2SpectBlue                          *SampledSpectrum
//)

type CoefficientSpectrum []float64

type RGBSpectrum struct {
	*CoefficientSpectrum
}
type SampledSpectrum struct {
	*CoefficientSpectrum
}

type DefaultSpectrum = RGBSpectrum

//type Spectrum = RGBSpectrum

type SpectrumGetter interface {
	Index(i int) float64
}

type SpectrumSetter interface {
	SetIndex(i int, v float64)
}

type Spectrum interface {
	SpectrumGetter
	SpectrumSetter

	Clone() Spectrum

	Add(other SpectrumGetter) Spectrum
	AddAssign(other SpectrumGetter) Spectrum
	AddScalar(other float64) Spectrum

	Sub(other SpectrumGetter) Spectrum
	SubAssign(other SpectrumGetter) Spectrum
	SubScalar(other float64) Spectrum

	Mul(other SpectrumGetter) Spectrum
	MulAssign(other SpectrumGetter) Spectrum
	MulScalar(other float64) Spectrum

	Div(other SpectrumGetter) Spectrum
	DivAssign(other SpectrumGetter) Spectrum
	DivScalar(other float64) Spectrum

	IsBlack() bool
	HasNaNs() bool
	Clamp(low, high float64) Spectrum
	Equals(other SpectrumGetter) bool

	SetAll(v float64) Spectrum

	Y() float64
}

func NewSpectrum(def float64) Spectrum {
	var s Spectrum = &DefaultSpectrum{}
	switch s.(type) {
	case RGBSpectrum:
		return NewRGBSpectrum(def, def, def)
	case SampledSpectrum:
		return NewSampledSpectrum(def)
	default:
		return NewRGBSpectrum(def, def, def)
	}
}

func NewRGBSpectrum(r, g, b float64) *RGBSpectrum {
	return &RGBSpectrum{&CoefficientSpectrum{r, g, b}}
}

func NewSampledSpectrum(d float64) *SampledSpectrum {
	cs := CoefficientSpectrum{}
	for i := 0; i < nSpectralSamples; i++ {
		cs = append(cs, d)
	}
	return &SampledSpectrum{&cs}
}

func (s CoefficientSpectrum) String() string {
	str := "["
	for i, v := range s {
		str += fmt.Sprintf("%v", v)
		if i < len(s) - 1 {
			str += ","
		}
	}
	str += "]"
	return str
}

func (s CoefficientSpectrum) Index(i int) float64 {
	return s[i]
}

func (s CoefficientSpectrum) SetIndex(i int, v float64) {
	s[i] = v
}

func (s CoefficientSpectrum) SetAll(v float64) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] = v
	}
	return s
}

func (s CoefficientSpectrum) Clone() Spectrum {
	ns := make(CoefficientSpectrum, len(s))
	copy(ns, s)
	return ns
}

func (s CoefficientSpectrum) Add(other SpectrumGetter) Spectrum {
	ns := s.Clone()
	ns.AddAssign(other)
	return ns
}

func (s CoefficientSpectrum) AddAssign(other SpectrumGetter) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] += other.Index(i)
	}
	return s
}

func (s CoefficientSpectrum) AddScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i) + other)
	}
	return ns
}

func (s CoefficientSpectrum) Sub(other SpectrumGetter) Spectrum {
	ns := s.Clone()
	ns.SubAssign(other)
	return ns
}

func (s CoefficientSpectrum) SubAssign(other SpectrumGetter) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] -= other.Index(i)
	}
	return s
}

func (s CoefficientSpectrum) SubScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i) - other)
	}
	return ns
}

func (s CoefficientSpectrum) Mul(other SpectrumGetter) Spectrum {
	ns := s.Clone()
	ns.MulAssign(other)
	return ns
}

func (s CoefficientSpectrum) MulAssign(other SpectrumGetter) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] *= other.Index(i)
	}
	return s
}

func (s CoefficientSpectrum) MulScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i) * other)
	}
	return ns
}

func (s CoefficientSpectrum) Div(other SpectrumGetter) Spectrum {
	ns := s.Clone()
	ns.DivAssign(other)
	return ns
}

func (s CoefficientSpectrum) DivAssign(other SpectrumGetter) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] /= other.Index(i)
	}
	return s
}

func (s CoefficientSpectrum) DivScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i) / other)
	}
	return ns
}

func (s CoefficientSpectrum) Clamp(low, high float64) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] = Clamp(s[i], low, high)
	}
	return s
}

func (s CoefficientSpectrum) Equals(other SpectrumGetter) bool {
	for i := 0; i < len(s); i++ {
		if s[i] != other.Index(i) {
			return false
		}
	}
	return true
}

func (s CoefficientSpectrum) IsBlack() bool {
	for i := 0; i < len(s); i++ {
		if s[i] != 0.0 {
			return false
		}
	}
	return true
}

func (s CoefficientSpectrum) HasNaNs() bool {
	for i := 0; i < len(s); i++ {
		if math.IsNaN(s[i]) {
			return true
		}
	}
	return false
}

func (s CoefficientSpectrum) Y() float64 {
	return 0
}

func (s RGBSpectrum) R() float64 {
	return s.Index(0)
}

func (s RGBSpectrum) G() float64 {
	return s.Index(1)
}

func (s RGBSpectrum) B() float64 {
	return s.Index(1)
}

func (s RGBSpectrum) Y() float64 {
	YWeight := [3]float64{0.212671, 0.715160, 0.072169}
	return YWeight[0]*s.R() + YWeight[1]*s.G() + YWeight[2]*s.B()
}