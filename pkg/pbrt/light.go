//go:generate mockgen -source=light.go -destination=light.mock.go -package=pbrt

package pbrt

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
	Preprocess(s Scene)

	SampleLi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, vis *VisibilityTester)
	Power() Spectrum
	Le(r *Ray) Spectrum
	PdfLi(ref Interaction, wi *Vector3f) float64
	SampleLe(u1, u2 *Point2f, time float64) (s Spectrum, r *Ray, nLight *Normal3f, pdfPos, pdfDir float64)
	PdfLe(r *Ray, nLight *Normal3f) (pdfPos, pdfDir float64)
}

func Le(r *Ray) Spectrum {
	return NewSpectrum(0.0)
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
	Flags                      LightFlag
	Samples                    int32
	MediumAccessor             *MediumAccessor
	LightToWorld, WorldToLight *Transform
}

func NewAreaLight(lightToWorld *Transform, mediumAccessor *MediumAccessor, nSamples int32) *AreaLight {
	return &AreaLight{
		Flags:          LightFlagArea,
		Samples:        nSamples,
		MediumAccessor: mediumAccessor,
		LightToWorld:   lightToWorld,
		WorldToLight:   lightToWorld.Inverse(),
	}
}

func (l *AreaLight) GetFlags() LightFlag {
	return l.Flags
}

func (l *AreaLight) GetSamples() int32 {
	return l.Samples
}

func (l *AreaLight) Preprocess(s Scene) {

}
