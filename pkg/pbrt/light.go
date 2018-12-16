//go:generate mockgen -source=light.go -destination=light.mock.go -package=pbrt

package pbrt

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
)

type LightFlag int

const (
	LightFlagDeltaPosition = 1 << iota
	LightFlagDeltaDirection
	LightFlagArea
	LightFlagInfinite
)

func IsDeltaLight(flags LightFlag) bool {
	return flags&LightFlagDeltaPosition >= 1 || flags&LightFlagDeltaDirection >= 1
}

type Light interface {
	GetFlags() LightFlag
	GetSamples() int32
	Preprocess()

	SampleLi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, vis *VisibilityTester)
	Power() Spectrum
	Le(r *RayDifferential) Spectrum
	PdfLi(ref Interaction, wi *Vector3f) float64
	SampleLe(u1, u2 *Point2f, time float64) (s Spectrum, r *Ray, nLight *Normal3f, pdfPos, pdfDir float64)
	PdfLe(r *Ray, nLight *Normal3f) (pdfPos, pdfDir float64)
}

func Le(r *RayDifferential) Spectrum {
	return NewSpectrum(0.0)
}

type PointLight struct {
	flags          LightFlag
	nSamples       int32
	MediumAccessor *MediumAccessor
	lightToWorld   *Transform
	worldToLight   *Transform

	pLight *Point3f
	I      Spectrum
}

func NewPointLight(lightToWorld *Transform, mediumAccessor *MediumAccessor, i Spectrum) *PointLight {
	p, _ := lightToWorld.TransformPoint(new(Point3f), new(Vector3f))
	return &PointLight{
		flags:          LightFlagDeltaPosition,
		nSamples:       8,
		MediumAccessor: mediumAccessor,
		lightToWorld:   lightToWorld,
		worldToLight:   lightToWorld.Inverse(),
		pLight:         p,
		I:              i,
	}
}

func (l *PointLight) GetFlags() LightFlag {
	return l.flags
}

func (l *PointLight) GetSamples() int32 {
	return l.nSamples
}

func (l *PointLight) Preprocess() {

}

func (l *PointLight) SampleLi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, vis *VisibilityTester) {
	wi = l.pLight.Sub(ref.GetPoint()).Normalized()
	pdf = 1.0
	vis = NewVisibilityTester(ref, NewInteraction(l.pLight, new(Vector3f), new(Vector3f), new(Vector3f), ref.GetTime(), l.MediumAccessor))
	return l.I.DivScalar(l.pLight.DistanceSquared(ref.GetPoint())), wi, pdf, vis
}

func (l *PointLight) Power() Spectrum {
	return l.I.MulScalar(4 * math.Pi)
}

func (l *PointLight) Le(r *RayDifferential) Spectrum {
	return Le(r)
}

func (l *PointLight) PdfLi(ref Interaction, wi *Vector3f) float64 {
	return 0
}

func (l *PointLight) SampleLe(u1, u2 *Point2f, time float64) (s Spectrum, r *Ray, nLight *Normal3f, pdfPos, pdfDir float64) {
	r = NewRayWithMedium(l.pLight, UniformSampleSphere(u1), time, l.MediumAccessor.Inside)
	return l.I, r, r.Direction, 1, UniformSpherePdf()
}

func (l *PointLight) PdfLe(r *Ray, nLight *Normal3f) (pdfPos, pdfDir float64) {
	return 0, UniformSpherePdf()
}

type VisibilityTester struct {
	p0, p1 Interaction
}

func NewVisibilityTester(p0, p1 Interaction) *VisibilityTester {
	return &VisibilityTester{
		p0: p0,
		p1: p1,
	}
}

func (v *VisibilityTester) Unoccluded(s Scene) bool {
	return !s.IntersectP(v.p0.SpawnRayToInteraction(v.p1))
}

func (v *VisibilityTester) Tr(scene Scene, sampler Sampler) Spectrum {
	ray := v.p0.SpawnRayToInteraction(v.p1)
	Tr := NewSpectrum(1)
	for {
		isect := NewSurfaceInteraction()
		hitSurface := scene.Intersect(ray, isect)
		// handle opaque surface along ray's path
		if hitSurface && isect.Primitive.GetMaterial() != nil {
			return NewSpectrum(0)
		}

		// update transmittance for current ray segment
		if ray.Medium != nil {
			Tr.MulAssign(ray.Medium.Tr(ray, sampler))
		}

		// generate next ray segment or return final transmittance
		if !hitSurface {
			break
		}
		ray = isect.SpawnRayToInteraction(v.p1)
	}
	return Tr
}

type AreaLighter interface {
	Light

	L(i Interaction, w *Vector3f) Spectrum
}

type AreaLight struct {
	flags                      LightFlag
	nSamples                   int
	MediumAccessor             *MediumAccessor
	lightToWorld, worldToLight *Transform
}

func NewAreaLight(lightToWorld *Transform, mediumAccessor *MediumAccessor, nSamples int) *AreaLight {
	return &AreaLight{
		flags:          LightFlagArea,
		nSamples:       nSamples,
		MediumAccessor: mediumAccessor,
		lightToWorld:   lightToWorld,
		worldToLight:   lightToWorld.Inverse(),
	}
}

// TODO: move DiffuseAreaLight out of core pbrt

type DiffuseAreaLight struct {
	*AreaLight

	LEmit    Spectrum
	shape    Shape
	twoSided bool
	area     float64
}

func NewDiffuseAreaLight(lightToWorld *Transform, mediumAccessor *MediumAccessor, LEmit Spectrum, nSamples int, shape Shape, twoSided bool) *DiffuseAreaLight {
	return &DiffuseAreaLight{
		AreaLight: NewAreaLight(lightToWorld, mediumAccessor, nSamples),
		LEmit:     LEmit,
		shape:     shape,
		twoSided:  twoSided,
		area:      shape.Area(),
	}
}

func (l *DiffuseAreaLight) Power() Spectrum {
	multiplier := 1.0
	if l.twoSided {
		multiplier = 2
	}
	return l.LEmit.MulScalar(multiplier * l.area * math.Pi)

}

func (l *DiffuseAreaLight) L(intr Interaction, w *Vector3f) Spectrum {
	if l.twoSided || intr.GetNormal().Dot(w) > 0 {
		return l.LEmit
	}
	return NewSpectrum(0)
}

func (l *DiffuseAreaLight) SampleLi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, vis *VisibilityTester) {
	pShape, pdf := l.shape.SampleAtInteraction(ref, u)
	pShape.SetMediumAccessor(l.MediumAccessor)
	if pdf == 0 || pShape.GetPoint().Sub(ref.GetPoint()).LengthSquared() == 0 {
		return NewSpectrum(0), new(Vector3f), 0, new(VisibilityTester)
	}

	wi = pShape.GetPoint().Sub(ref.GetPoint())
	vis = NewVisibilityTester(ref, pShape)

	return l.L(pShape, wi.MulScalar(-1)), wi, pdf, vis
}

func (l *DiffuseAreaLight) PdfLi(ref Interaction, wi *Vector3f) float64 {
	return l.shape.PdfWi(ref, wi)
}

func (l *DiffuseAreaLight) SampleLe(u1, u2 *Point2f, time float64) (s Spectrum, r *Ray, nLight *Normal3f, pdfPos, pdfDir float64) {
	// sample a Point on the are light's shape, pShape
	pShape, pdfPos := l.shape.Sample(u1)
	pShape.SetMediumAccessor(l.MediumAccessor)

	var w *Vector3f
	if l.twoSided {
		u := u2
		// choose a side to sample and then remap u[0] to [0,1] before
		// applying cosine-weighted hemispher sampling for the chose side.
		if u.X < 0.5 {
			u.X = math.Min(u.X*2, math.OneMinusEpsilon)
			w = CosineSampleHemisphere(u)
		} else {
			u.X = math.Min((u.X-0.5)*2, math.OneMinusEpsilon)
			w = CosineSampleHemisphere(u)
			w.Z *= -1
		}
		pdfDir = 0.5 * CosineHemispherePdf(math.Abs(w.Z))
	} else {
		w = CosineSampleHemisphere(u2)
		pdfDir = CosineHemispherePdf(w.Z)
	}

	v1, v2 := CoordinateSystem(pShape.GetNormal())
	w = v1.MulScalar(w.X).Add(v2.MulScalar(w.Y)).Add(pShape.GetNormal().MulScalar(w.Z))
	return l.L(pShape, w), pShape.SpawnRay(w), pShape.GetNormal(), pdfPos, pdfDir
}

func (l *DiffuseAreaLight) PdfLe(r *Ray, n *Normal3f) (pdfPos, pdfDir float64) {
	it := NewInteraction(r.Origin, &Vector3f{}, n, n, r.Time, l.MediumAccessor)
	pdfPos = l.shape.Pdf(it)
	if l.twoSided {
		pdfPos = 0.5 * CosineHemispherePdf(n.AbsDot(r.Direction))
	} else {
		pdfPos = CosineHemispherePdf(n.AbsDot(r.Direction))
	}
	return pdfPos, pdfDir
}
