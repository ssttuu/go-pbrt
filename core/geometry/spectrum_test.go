package geometry

import (
	"testing"
)

func TestSpectrumMultiply(t *testing.T) {
	s1 := NewRGBSpectrum([]float64{1, 2, 3})
	s2 := NewRGBSpectrum([]float64{2, 3, 4})

	s3 := SpectrumAdd(&s1, &s2)
	expected := NewRGBSpectrum([]float64{3, 5, 7})
	if SpectrumNotEquals(s3, &expected) {
		t.Fail()
	}
}