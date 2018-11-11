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

func NewSphere(name string, objectToWorld, worldToObject *Transform, reverseOrientation bool, radius, zMin, zMax, phiMax float64) *Sphere {
	return &Sphere{
		Shape: &Shape{
			name:               name,
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

func NewSphereShape(name string, o2w, w2o *Transform, reverseOrientation bool, radius float64) Shaper {
	return NewSphere(name, o2w, w2o, reverseOrientation, radius, -radius, radius, 360.0)
}

func (s *Sphere) ObjectBound() *Bounds3 {
	return &Bounds3{
		Min: &Point3f{-s.radius, -s.radius, s.zMin},
		Max: &Point3f{s.radius, s.radius, s.zMax},
	}
}

func (s *Sphere) WorldBound() *Bounds3 {
	return s.objectToWorld.TransformBounds(s.ObjectBound())
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

func (s *Sphere) Sample(u *Point2f) (i Interaction, pdf float64) {
	pObj := new(Point3f).Add(UniformSampleSphere(u).MulScalar(s.radius))

	var it *interaction
	it.normal = s.objectToWorld.TransformNormal((*Normal3f)(pObj)).Normalized()
	if s.reverseOrientation {
		it.normal = it.normal.MulScalar(-1)
	}

	// reproject pObj to sphere surface
	pObj = pObj.MulScalar(s.radius / pObj.Distance(&Point3f{0, 0, 0}))
	it.point = s.objectToWorld.TransformPoint(pObj)

	return it, 1.0 / s.Area()
}

func (s *Sphere) SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64) {
	pCenter := s.objectToWorld.TransformPoint(&Point3f{0, 0, 0})

	pOrigin := ref.GetPoint()
	if pOrigin.Distance(pCenter) <= s.radius {
		intr, pdf := s.Sample(u)
		wi := intr.GetPoint().Sub(ref.GetPoint())
		if wi.LengthSquared() == 0 {
			pdf = 0
		} else {
			wi.Normalize()
			pdf = ref.GetPoint().LengthSquared() / intr.GetNormal().AbsDot(wi.MulScalar(-1))
		}

		if math.IsInf(pdf, 1) {
			pdf = 0.0
		}
		return intr, pdf
	}

	// compute coordinate system for sphere sampling
	wc := pCenter.Sub(ref.GetPoint()).Normalized()
	wcX, wcY := CoordinateSystem(wc)

	radiusSquared := s.radius * s.radius

	// sample sphere uniformly inside subtended code

	// compute theta and phi values for sample in code
	sinThetaMax2 := radiusSquared / ref.GetPoint().DistanceSquared(pCenter)
	cosThetaMax := math.Sqrt(math.Max(0, 1.0-sinThetaMax2))
	cosTheta := (1.0 - u.X) + u.X*cosThetaMax
	sinTheta := math.Sqrt(math.Max(0, 1-cosTheta*cosTheta))
	phi := u.Y * 2 * math.Pi

	// compute angle alpha from center of sphere to sampled point on surface
	dc := ref.GetPoint().Distance(pCenter)
	ds := dc*cosTheta - math.Sqrt(math.Max(0, radiusSquared-dc*dc*sinTheta*sinTheta))
	cosAlpha := (dc*dc + radiusSquared - ds*ds) / (2.0 * dc * s.radius)
	sinAlpha := math.Sqrt(math.Max(0, 1.0-cosAlpha*cosAlpha))

	// compute surface normal and sampled point on sphere
	nWorld := SphericalDirectionXYZ(sinAlpha, cosAlpha, phi, wcX.MulScalar(-1), wcY.MulScalar(-1), wc.MulScalar(-1))
	pWorld := pCenter.Add(nWorld.MulScalar(s.radius))

	// return interaction for sampled point on sphere
	it := new(interaction)
	it.point = pWorld
	it.normal = (*Normal3f)(nWorld)
	if s.reverseOrientation {
		it.normal = it.normal.MulScalar(-1)
	}

	// uniform code PDF
	pdf = 1.0 / (2.0 * Pi * (1 - cosThetaMax))

	return it, pdf
}

func (s *Sphere) PdfWi(ref Interaction) (float64, *Vector3f) {
	pCenter := s.objectToWorld.TransformPoint(&Point3f{0, 0, 0})
	pOrigin := ref.GetPoint()

	if pOrigin.Distance(pCenter) <= s.radius {
		return s.Shape.PdfWi(ref)
	}

	sinThetaMax2 := s.radius * s.radius / ref.GetPoint().Dot(pCenter)
	cosThetaMax := math.Sqrt(math.Max(0, 1.0-sinThetaMax2))

	return UniformConePdf(cosThetaMax), new(Vector3f)
}

func (s *Sphere) SolidAngle(p *Point3f, nSamples int) float64 {
	pCenter := s.objectToWorld.TransformPoint(&Point3f{0, 0, 0})
	if p.Distance(pCenter) <= s.radius {
		return 4 * math.Pi
	}

	sinTheta2 := s.radius * s.radius / p.Dot(pCenter)
	cosTheta := math.Sqrt(math.Max(0, 1.0-sinTheta2))

	return 2 * math.Pi * (1.0 - cosTheta)
}
