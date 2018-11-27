package pbrt

import (
	"github.com/stupschwartz/go-pbrt/pkg/efloat"
	"github.com/stupschwartz/go-pbrt/pkg/math"
)

type Sphere struct {
	name string

	objectToWorld, worldToObject *Transform
	reverseOrientation           bool
	transformSwapsHandedness     bool

	radius                     float64
	zMin, zMax                 float64
	thetaMin, thetaMax, phiMax float64
}

func NewSphere(name string, objectToWorld, worldToObject *Transform, reverseOrientation bool, radius, zMin, zMax, phiMax float64) *Sphere {
	return &Sphere{
		name:               name,
		objectToWorld:      objectToWorld,
		worldToObject:      worldToObject,
		reverseOrientation: reverseOrientation,
		radius:             radius,
		zMin:               math.Clamp(math.Min(zMin, zMax), -radius, radius),
		zMax:               math.Clamp(math.Max(zMin, zMax), -radius, radius),
		thetaMin:           math.Acos(math.Clamp(math.Min(zMin, zMax)/radius, -1, 1)),
		thetaMax:           math.Acos(math.Clamp(math.Max(zMin, zMax)/radius, -1, 1)),
		phiMax:             math.Radians(math.Clamp(phiMax, 0, 360)),
	}
}

func NewSphereShape(name string, o2w *Transform, reverseOrientation bool, radius float64) Shape {
	return NewSphere(name, o2w, o2w.Inverse(), reverseOrientation, radius, -radius, radius, 360.0)
}

func (s *Sphere) GetName() string {
	return s.name
}

func (s *Sphere) Area() float64 {
	return s.phiMax * s.radius * (s.zMax - s.zMin)
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

func (s *Sphere) ReverseOrientation() bool {
	return s.reverseOrientation
}
func (s *Sphere) TransformSwapsHandedness() bool {
	return s.transformSwapsHandedness
}

func (s *Sphere) Intersect(r *Ray, si *SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64) {
	ray, originError, directionError := s.worldToObject.TransformRay(r)

	// Initialize _EFloat_ ray coordinate values
	ox := efloat.New(ray.Origin.X, originError.X)
	oy := efloat.New(ray.Origin.Y, originError.Y)
	oz := efloat.New(ray.Origin.Z, originError.Z)

	dx := efloat.New(ray.Direction.X, directionError.X)
	dy := efloat.New(ray.Direction.Y, directionError.Y)
	dz := efloat.New(ray.Direction.Z, directionError.Z)

	a := dx.Mul(dx).Add(dy.Mul(dy)).Add(dz.Mul(dz))
	b := dx.Mul(ox).Add(dy.Mul(oy)).Add(dz.Mul(oz)).MulScalar(2.0)
	c := ox.Mul(ox).Add(oy.Mul(oy)).Add(oz.Mul(oz)).Sub(efloat.New(s.radius, 0).MulScalar(s.radius))

	t0, t1, ok := efloat.Quadratic(a, b, c)
	if !ok {
		return false, 0
	}

	// Check quadratic shape t0 and t1 for nearest inersection
	if t0.High > ray.TMax || t1.Low <= 0 {
		return false, 0
	}
	tShapeHit := t0
	if tShapeHit.Low <= 0 {
		tShapeHit = t1
		if tShapeHit.High > ray.TMax {
			return false, 0
		}
	}

	// compute sphere hit position
	pHit := ray.PointAt(tShapeHit.Value)

	// refine sphere intersection point
	pHit = pHit.MulScalar(s.radius / pHit.Distance(new(Point3f)))
	if pHit.X == 0.0 && pHit.Y == 0.0 {
		pHit.X = 1e-5 * s.radius
	}
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0.0 {
		phi += 2 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
		if tShapeHit == t1 {
			return false, 0
		}
		if t1.High > ray.TMax {
			return false, 0
		}

		tShapeHit = t1
		pHit = ray.PointAt(tShapeHit.Value)

		// refine sphere intersection Point
		pHit = pHit.MulScalar(s.radius / pHit.Distance(new(Point3f)))
		if pHit.X == 0.0 && pHit.Y == 0.0 {
			pHit.X = 1e-5 * s.radius
		}
		phi := math.Atan2(pHit.Y, pHit.X)
		if phi < 0.0 {
			phi += 2 * math.Pi
		}

		if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
			return false, 0
		}
	}

	// Find parametric representation of sphere hit
	u := phi / s.phiMax
	theta := math.Acos(math.Clamp(pHit.Z/s.radius, -1, 1))
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
	d2Pdvv := pHit.MulScalar(-(s.thetaMax - s.thetaMin) * (s.thetaMax - s.thetaMin))

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
	dndu := dpdu.MulScalar((f*F - e*G) * invEGF2).Add(dpdv.MulScalar((e*F - f*E) * invEGF2))
	dndv := dpdu.MulScalar((g*F - f*G) * invEGF2).Add(dpdv.MulScalar((f*F - g*E) * invEGF2))

	// compute error bounds for sphere inertsection
	pError := pHit.Abs().MulScalar(math.Gamma(5))

	*si = *NewSurfaceInteractionWith(pHit, pError, &Point2f{u, v}, ray.Direction.MulScalar(-1), dpdu, dpdv, dndu, dndv, ray.Time, s, 0)
	*si = *s.objectToWorld.TransformSurfaceInteraction(si)

	return true, tShapeHit.Value
}

func (s *Sphere) IntersectP(r *Ray, testAlphaTexture bool) (intersects bool) {
	// Transform _Ray_ to object space
	ray, originError, directionError := s.worldToObject.TransformRay(r)

	// Compute quadratic sphere coefficients

	// Initialize _EFloat_ ray coordinate values
	ox := efloat.New(ray.Origin.X, originError.X)
	oy := efloat.New(ray.Origin.Y, originError.Y)
	oz := efloat.New(ray.Origin.Z, originError.Z)

	dx := efloat.New(ray.Direction.X, directionError.X)
	dy := efloat.New(ray.Direction.Y, directionError.Y)
	dz := efloat.New(ray.Direction.Z, directionError.Z)

	a := dx.Mul(dx).Add(dy.Mul(dy)).Add(dz.Mul(dz))
	b := dx.Mul(ox).Add(dy.Mul(oy)).Add(dz.Mul(oz)).MulScalar(2.0)
	c := ox.Mul(ox).Add(oy.Mul(oy)).Add(oz.Mul(oz)).Sub(efloat.New(s.radius, 0).MulScalar(s.radius))

	// Solve quadratic equation for _t_ values
	t0, t1, ok := efloat.Quadratic(a, b, c)
	if !ok {
		return false
	}

	if t0.High > ray.TMax || t1.Low <= 0 {
		return false
	}
	tShapeHit := t0
	if tShapeHit.Low <= 0 {
		tShapeHit = t1
		if tShapeHit.High > ray.TMax {
			return false
		}
	}

	// compute sphere hit position and phi
	pHit := ray.PointAt(tShapeHit.Value)

	// refine sphere intersection Point
	pHit = pHit.MulScalar(s.radius / pHit.Distance(&Point3f{0, 0, 0}))
	if pHit.X == 0.0 && pHit.Y == 0.0 {
		pHit.X = 1e-5 * s.radius
	}
	phi := math.Atan2(pHit.Y, pHit.X)
	if phi < 0.0 {
		phi += 2 * math.Pi
	}

	// Test sphere intersection against clipping parameters
	if (s.zMin > -s.radius && pHit.Z < s.zMin) || s.zMax < s.radius && pHit.Z > s.zMax || phi > s.phiMax {
		if tShapeHit == t1 {
			return false
		}
		if t1.High > ray.TMax {
			return false
		}

		// compute sphere hit position and phi
		tShapeHit = t1
		pHit = ray.PointAt(tShapeHit.Value)

		// refine sphere intersection Point
		pHit = pHit.MulScalar(s.radius / pHit.Distance(new(Point3f)))
		if pHit.X == 0.0 && pHit.Y == 0.0 {
			pHit.X = 1e-5 * s.radius
		}
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
	pObj := UniformSampleSphere(u).MulScalar(s.radius)

	var it *interaction
	it.Normal = s.objectToWorld.TransformNormal(pObj).Normalized()
	if s.reverseOrientation {
		it.Normal = it.Normal.MulScalar(-1)
	}

	// reproject pObj to sphere surface
	pObj = pObj.MulScalar(s.radius / pObj.Distance(new(Point3f)))
	pObjError := pObj.Abs().MulScalar(math.Gamma(5))
	it.Point, it.PointError = s.objectToWorld.TransformPoint(pObj, pObjError)

	return it, 1.0 / s.Area()
}

func (s *Sphere) SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64) {
	pCenter, _ := s.objectToWorld.TransformPoint(new(Point3f), new(Vector3f))

	// Sample uniformly on sphere if $\pt{}$ is inside it
	pOrigin := OffsetRayOrigin(ref.GetPoint(), ref.GetPointError(), ref.GetNormal(), pCenter.Sub(ref.GetPoint()))
	if pOrigin.DistanceSquared(pCenter) <= s.radius*s.radius {
		intr, pdf := s.Sample(u)
		wi := intr.GetPoint().Sub(ref.GetPoint())
		if wi.LengthSquared() == 0 {
			pdf = 0
		} else {
			// Convert from area measure returned by Sample() call above to solid angle measure.
			wi.Normalize()
			pdf *= ref.GetPoint().DistanceSquared(intr.GetPoint()) / intr.GetNormal().AbsDot(wi.MulScalar(-1))
		}

		if math.IsInf(pdf) {
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

	// compute angle alpha from center of sphere to sampled Point on surface
	dc := ref.GetPoint().Distance(pCenter)
	ds := dc*cosTheta - math.Sqrt(math.Max(0, radiusSquared-(dc*dc)*(sinTheta*sinTheta)))
	cosAlpha := (dc*dc + radiusSquared - ds*ds) / (2.0 * dc * s.radius)
	sinAlpha := math.Sqrt(math.Max(0, 1.0-cosAlpha*cosAlpha))

	// compute surface Normal and sampled Point on sphere
	nWorld := SphericalDirectionXYZ(sinAlpha, cosAlpha, phi, wcX.MulScalar(-1), wcY.MulScalar(-1), wc.MulScalar(-1))
	pWorld := pCenter.Add(nWorld.MulScalar(s.radius))

	// return interaction for sampled Point on sphere
	it := new(interaction)
	it.Point = pWorld
	it.PointError = pWorld.Abs().MulScalar(math.Gamma(5.0))
	it.Normal = nWorld
	if s.reverseOrientation {
		it.Normal = it.Normal.MulScalar(-1)
	}

	return it, UniformConePdf(cosThetaMax)
}

func (s *Sphere) Pdf(ref Interaction) float64 {
	return Pdf(s, ref)
}

func (s *Sphere) PdfWi(ref Interaction, wi *Vector3f) float64 {
	pCenter, _ := s.objectToWorld.TransformPoint(new(Point3f), new(Vector3f))
	pOrigin := OffsetRayOrigin(ref.GetPoint(), ref.GetPointError(), ref.GetNormal(), ref.GetPoint().Sub(pCenter))

	if pOrigin.DistanceSquared(pCenter) <= s.radius*s.radius {
		return PdfWi(s, ref, wi)
	}

	// Compute general sphere PDF
	sinThetaMax2 := s.radius * s.radius / ref.GetPoint().DistanceSquared(pCenter)
	cosThetaMax := math.Sqrt(math.Max(0, 1.0-sinThetaMax2))

	return UniformConePdf(cosThetaMax)
}

func (s *Sphere) SolidAngle(p *Point3f, nSamples int) float64 {
	pCenter, _ := s.objectToWorld.TransformPoint(new(Point3f), new(Vector3f))
	if p.DistanceSquared(pCenter) <= s.radius*s.radius {
		return 4 * math.Pi
	}

	sinTheta2 := s.radius * s.radius / p.DistanceSquared(pCenter)
	cosTheta := math.Sqrt(math.Max(0, 1.0-sinTheta2))

	return 2 * math.Pi * (1.0 - cosTheta)
}
