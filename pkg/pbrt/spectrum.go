package pbrt

type Spectrum []float64

func NewSpectrum(samples ...float64) Spectrum {
	s := append(Spectrum{}, samples...)
	return s
}

func NewRGBSpectrum(r, g, b float64) Spectrum {
	return Spectrum{r, g, b}
}

func (s Spectrum) Add(other Spectrum) {
	for i := 0; i < len(s); i++ {
		s[i] += other[i]
	}
}

func (s Spectrum) Sub(other Spectrum) {
	for i := 0; i < len(s); i++ {
		s[i] -= other[i]
	}
}

func (s Spectrum) Mul(other Spectrum) {
	for i := 0; i < len(s); i++ {
		s[i] *= other[i]
	}
}

func (s Spectrum) Div(other Spectrum) {
	for i := 0; i < len(s); i++ {
		s[i] /= other[i]
	}
}

func (s Spectrum) Equals(other Spectrum) bool {
	for i := 0; i < len(s); i++ {
		if s[i] != other[i] {
			return false
		}
	}
	return true
}

