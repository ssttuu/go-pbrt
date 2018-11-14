package pbrt

import (
	"math"
	"fmt"
)

const (
	nSpectralSamples  = 60
	nRGB2SpectSamples = 32
)

type SpectrumType string

const (
	SpectrumTypeRGB     SpectrumType = "RGB"
	SpectrumTypeSampled              = "Sampled"
)

const DefaultSpectrumType = SpectrumTypeRGB

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

type Spectrum []float64

func NewSpectrum(def float64) Spectrum {
	switch DefaultSpectrumType {
	case SpectrumTypeRGB:
		return NewRGBSpectrum(def, def, def)
	case SpectrumTypeSampled:
		return NewSampledSpectrum(def)
	default:
		return NewRGBSpectrum(def, def, def)
	}
}

func NewRGBSpectrum(r, g, b float64) Spectrum {
	return Spectrum{r, g, b}
}

func NewSampledSpectrum(d float64) Spectrum {
	cs := Spectrum{}
	for i := 0; i < nSpectralSamples; i++ {
		cs = append(cs, d)
	}
	return cs
}

func (s Spectrum) String() string {
	str := "["
	for i, v := range s {
		str += fmt.Sprintf("%v", v)
		if i < len(s)-1 {
			str += ","
		}
	}
	str += "]"
	return str
}

func (s Spectrum) Index(i int) float64 {
	return s[i]
}

func (s Spectrum) SetIndex(i int, v float64) {
	s[i] = v
}

func (s Spectrum) SetAll(v float64) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] = v
	}
	return s
}

func (s Spectrum) Clone() Spectrum {
	ns := make(Spectrum, len(s))
	copy(ns, s)
	return ns
}

func (s Spectrum) Add(other Spectrum) Spectrum {
	ns := s.Clone()
	ns.AddAssign(other)
	return ns
}

func (s Spectrum) AddAssign(other Spectrum) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] += other.Index(i)
	}
	return s
}

func (s Spectrum) AddScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i)+other)
	}
	return ns
}

func (s Spectrum) Sub(other Spectrum) Spectrum {
	ns := s.Clone()
	ns.SubAssign(other)
	return ns
}

func (s Spectrum) SubAssign(other Spectrum) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] -= other.Index(i)
	}
	return s
}

func (s Spectrum) SubScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i)-other)
	}
	return ns
}

func (s Spectrum) Mul(other Spectrum) Spectrum {
	ns := s.Clone()
	ns.MulAssign(other)
	return ns
}

func (s Spectrum) MulAssign(other Spectrum) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] *= other.Index(i)
	}
	return s
}

func (s Spectrum) MulScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i)*other)
	}
	return ns
}

func (s Spectrum) Div(other Spectrum) Spectrum {
	ns := s.Clone()
	ns.DivAssign(other)
	return ns
}

func (s Spectrum) DivAssign(other Spectrum) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] /= other.Index(i)
	}
	return s
}

func (s Spectrum) DivScalar(other float64) Spectrum {
	ns := s.Clone()
	for i := 0; i < len(s); i++ {
		ns.SetIndex(i, ns.Index(i)/other)
	}
	return ns
}

func (s Spectrum) Clamp(low, high float64) Spectrum {
	for i := 0; i < len(s); i++ {
		s[i] = Clamp(s[i], low, high)
	}
	return s
}

func (s Spectrum) Equals(other Spectrum) bool {
	for i := 0; i < len(s); i++ {
		if s[i] != other.Index(i) {
			return false
		}
	}
	return true
}

func (s Spectrum) IsBlack() bool {
	for i := 0; i < len(s); i++ {
		if s[i] != 0.0 {
			return false
		}
	}
	return true
}

func (s Spectrum) HasNaNs() bool {
	for i := 0; i < len(s); i++ {
		if math.IsNaN(s[i]) {
			return true
		}
	}
	return false
}

func (s Spectrum) Y() float64 {
	return 0
}
