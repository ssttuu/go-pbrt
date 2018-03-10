package pbrt

import (
	"testing"
)

func TestSpectrumMultiply(t *testing.T) {
	s1 := NewSpectrum(1)
	s2 := NewSpectrum(2)

	s1.AddAssign(s2)
	expected := NewSpectrum(3)
	if !s1.Equals(expected) {
		t.Errorf("%+v %+v %+v", s1, s2, expected)
	}
}
