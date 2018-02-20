package geometry

import "math"

type Sphere struct {
	Center *VectorXYZ
	Radius float64
	//RadiusSquared float64

	SurfaceColor *VectorXYZ
	EmissionColor *VectorXYZ

	Transparency float64
	Reflection float64
}

const INFINITY float64 = 10e8

func (s *Sphere) Intersect(origin, direction *VectorXYZ) (t0 float64, t1 float64, intersects bool) {
	intersects = false
	t0 = INFINITY
	t1 = INFINITY

	l := VectorXYZSubtract(s.Center, origin)
	tca := l.Dot(direction)
	if tca < 0.0 {
		intersects = false
		return
	}

	d2 := l.Dot(l) - tca * tca
	if d2 > (s.Radius * s.Radius) {
		intersects = false
		return
	}

	thc := math.Sqrt((s.Radius * s.Radius) - d2)

	t0 = tca - thc
	t1 = tca + thc
	intersects = true

	return
}