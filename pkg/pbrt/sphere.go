package pbrt

import (
	"math"
)

type Sphere struct {
	Center        *Vector3
	Radius        float64
	RadiusSquared float64

	SurfaceColor  *Vector3
	EmissionColor *Vector3

	Transparency  float64
	Reflection    float64
}

//func NewSphere(center *VectorXYZ, Radius float64)



func (s *Sphere) Intersects(origin, direction *Vector3) (float64, float64, bool) {
	l := s.Center.Sub(origin)
	tca := l.Dot(direction)

	if tca < 0.0 {
		return Infinity, Infinity, false
	}

	d2 := l.LengthSquared() - tca * tca
	if d2 > (s.Radius * s.Radius) {
		return Infinity, Infinity, false
	}

	thc := math.Sqrt(s.Radius * s.Radius) - d2

	return tca - thc, tca + thc, true
}