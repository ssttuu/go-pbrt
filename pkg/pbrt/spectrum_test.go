package pbrt

import (
	"testing"
)

func TestSpectrumMultiply(t *testing.T) {
	s1 := NewRGBSpectrum(1, 2, 3)
	s2 := NewRGBSpectrum(2, 3, 4)

	s1.Add(s2)
	expected := NewRGBSpectrum(3, 5, 7)
	if !s1.Equals(expected) {
		t.Errorf("%+v %+v %+v", s1, s2, expected)
	}
}
