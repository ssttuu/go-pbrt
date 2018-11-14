package pbrt

import (
	"testing"
	"github.com/stretchr/testify/assert"
)


func TestSpectrum_Add(t *testing.T) {
	s1 := NewSpectrum(1)
	s2 := NewSpectrum(2)

	result := s1.Add(s2)

	assert.Equal(t, NewSpectrum(1), s1)
	assert.Equal(t, NewSpectrum(3), result)
}

func TestSpectrum_AddAssign(t *testing.T) {
	s1 := NewSpectrum(1)
	s2 := NewSpectrum(2)

	result := s1.AddAssign(s2)

	assert.Equal(t, NewSpectrum(3), s1)
	assert.Equal(t, s1, result)
}

func TestSpectrum_AddScalar(t *testing.T) {
	s1 := NewSpectrum(1)

	result := s1.AddScalar(2.0)

	assert.Equal(t, NewSpectrum(1), s1)
	assert.Equal(t, NewSpectrum(3), result)
}

func TestSpectrum_IsBlack(t *testing.T) {
	assert.False(t, NewRGBSpectrum(1.0, 0.0, 1.0).IsBlack())
	assert.False(t, NewSpectrum(1).IsBlack())
	assert.False(t, NewSpectrum(0.0001).IsBlack())
	assert.True(t, NewSpectrum(0.0).IsBlack())
}