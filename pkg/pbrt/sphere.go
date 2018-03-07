package pbrt

import (
	"math"
)

type Sphere struct {
	*Shape

	radius                     float64
	zMin, zMax                 float64
	thetaMin, thetaMax, phiMax float64
}

func NewSphere(objectToWorld, worldToObject *Transform, reverseOrientation bool, radius, zMin, zMax, phiMax float64) *Sphere {
	return &Sphere{
		Shape: &Shape{
			objectToWorld:      objectToWorld,
			worldToObject:      worldToObject,
			reverseOrientation: reverseOrientation,
		},
		radius:   radius,
		zMin:     Clamp(math.Min(zMin, zMax), -radius, radius),
		zMax:     Clamp(math.Max(zMin, zMax), -radius, radius),
		thetaMin: math.Acos(Clamp(math.Min(zMin, zMax)/radius, -1, 1)),
		thetaMax: math.Acos(Clamp(math.Max(zMin, zMax)/radius, -1, 1)),
		phiMax:   Radians(Clamp(phiMax, 0, 360)),
	}
}

func NewSphereShape(o2w, w2o *Transform, reverseOrientation bool, radius float64) Shaper {
	return NewSphere(o2w, w2o, reverseOrientation, radius, -radius, radius, 360.0)
}

func (s *Sphere) ObjectBound() *Bounds3 {
	return &Bounds3{
		min: &Point3f{-s.radius, -s.radius, s.zMin},
		max: &Point3f{s.radius, s.radius, s.zMax},
	}
}

func (s *Sphere) Intersect(r *Ray, testAlphaTexture bool) (intersects bool, tHit float64, si *SurfaceInteraction) {
	ray := s.worldToObject.TransformRay(r)

	a := ray.direction.LengthSquared()
	b := 2 * ray.direction.Dot(ray.origin)
	c := ray.origin.LengthSquared() - s.radius*s.radius

	hasRoot, t0, t1 := Quadratic(a, b, c)
	if !hasRoot {
		return false, 0, nil
	}

	// skipped floating point error handling
	tShapeHit := t0

	// compute sphere hit position
	pHit := ray.PointAt(tShapeHit)

	pHit = pHit.MulScalar(s.radius / pHit.Distance(&Point3f{0, 0, 0}))
	if pHit.X == 0.0 && pHit.Y == 0.0 {
		pHit.X = 1e-5 * s.radius
	}

	// refine sphere intersection point
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0.0 {
		phi += 2 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
		tShapeHit = t1
		pHit = ray.PointAt(tShapeHit)

		// refine sphere intersection point
		pHit = pHit.MulScalar(s.radius / pHit.Distance(&Point3f{0, 0, 0}))
		if pHit.X == 0.0 && pHit.Y == 0.0 {
			pHit.X = 1e-5 * s.radius
		}

		// refine sphere intersection point
		phi := math.Atan2(pHit.Y, pHit.X)
		if phi < 0.0 {
			phi += 2 * math.Pi
		}

		if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
			return false, 0, nil
		}
	}

	// Find parametric representation of sphere hit
	u := phi / s.phiMax
	theta := math.Acos(Clamp(pHit.Z/s.radius, -1, 1))
	v := (theta - s.thetaMin) / (s.thetaMax - s.thetaMin)

	// compute sphere dpdu and dpdv
	zRadius := math.Sqrt(pHit.X*pHit.X + pHit.Y*pHit.Y)
	invZRadius := 1.0 / zRadius
	cosPhi := pHit.X * invZRadius
	sinPhi := pHit.Y * invZRadius
	dpdu := &Vector3f{-s.phiMax * pHit.Y, s.phiMax * pHit.X, 0}
	dpdv := new(Vector3f).Set(pHit.Z*cosPhi, pHit.Z*sinPhi, -s.radius*math.Sin(theta)).MulScalar(s.thetaMax - s.thetaMin)

	// compute sphere dndu and dndv
	d2Pduu := new(Vector3f).Set(pHit.X, pHit.Y, 0.0).MulScalar(-s.phiMax * s.phiMax)
	d2Pduv := new(Vector3f).Set(-sinPhi, cosPhi, 0.0).MulScalar((s.thetaMax - s.thetaMin) * pHit.Z * s.phiMax)
	d2Pdvv := new(Vector3f).Set(pHit.X, pHit.Y, pHit.Z).MulScalar(-(s.thetaMax - s.thetaMin) * (s.thetaMax - s.thetaMin))

	// compute coefficients for fundamental forms
	E := dpdu.Dot(dpdu)
	F := dpdu.Dot(dpdv)
	G := dpdv.Dot(dpdv)
	N := dpdu.Cross(dpdv).Normalized()
	e := N.Dot(d2Pduu)
	f := N.Dot(d2Pduv)
	g := N.Dot(d2Pdvv)

	//compute dndu and dndv from fundamental form coefficients
	invEGF2 := 1.0 / (E*G - F*F)
	dndu := (*Normal3f)(dpdu.MulScalar((f*F - e*G) * invEGF2).Add(dpdv.MulScalar((e*F - f*E) * invEGF2)))
	dndv := (*Normal3f)(dpdu.MulScalar((g*F - f*G) * invEGF2).Add(dpdv.MulScalar((f*F - g*E) * invEGF2)))

	// skipping error bounds for sphere intersection

	si = NewSurfaceInteraction(pHit, &Vector3f{}, &Point2f{u, v}, ray.direction.MulScalar(-1), dpdu, dpdv, dndu, dndv, ray.time, s, 0)
	si = s.objectToWorld.TransformSurfaceInteraction(si)

	return true, tShapeHit, si
}

func (s *Sphere) IntersectP(r *Ray, testAlphaTexture bool) (intersects bool) {
	ray := s.worldToObject.TransformRay(r)

	a := ray.direction.LengthSquared()
	b := 2 * ray.direction.Dot(ray.origin)
	c := ray.origin.LengthSquared() - s.radius*s.radius

	hasRoot, t0, t1 := Quadratic(a, b, c)
	if !hasRoot {
		return false
	}

	// skipped floating point error handling
	tShapeHit := t0

	// compute sphere hit position
	pHit := ray.PointAt(tShapeHit)

	pHit = pHit.MulScalar(s.radius / pHit.Distance(&Point3f{0, 0, 0}))
	if pHit.X == 0.0 && pHit.Y == 0.0 {
		pHit.X = 1e-5 * s.radius
	}

	// refine sphere intersection point
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0.0 {
		phi += 2 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
		tShapeHit = t1
		pHit = ray.PointAt(tShapeHit)

		// refine sphere intersection point
		pHit = pHit.MulScalar(s.radius / pHit.Distance(&Point3f{0, 0, 0}))
		if pHit.X == 0.0 && pHit.Y == 0.0 {
			pHit.X = 1e-5 * s.radius
		}

		// refine sphere intersection point
		phi := math.Atan2(pHit.Y, pHit.X)
		if phi < 0.0 {
			phi += 2 * math.Pi
		}

		if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
			return false
		}
	}

	return true
}
