package pbrt

import (
	"math"
)

type Interactioner interface {
	//IsSurfaceInteraction() bool
	//IsMediumInteraction() bool
	//
	SpawnRay(direction *Vector3f) *Ray
	GetPoint() *Point3f
	//GetPointError() *Vector3f
	//GetTime() float64
	GetNormal() *Normal3f
	//GetMedium(*Vector3f) *Mediumer
}

type Interaction struct {
	point          *Point3f
	pError         *Vector3f
	time           float64
	wo             *Vector3f
	normal         *Normal3f
	mediumAccessor *MediumAccessor
}

func (i *Interaction) GetPoint() *Point3f {
	return i.point
}

func (i *Interaction) GetNormal() *Normal3f {
	return i.normal
}

func (i *Interaction) SpawnRay(direction *Vector3f) *Ray {
	origin := OffsetRayOrigin(i.point, i.pError, i.normal, direction)
	return &Ray{origin, direction, Infinity, i.time, i.GetMedium(direction)}
}

func (i *Interaction) SpawnRayToPoint(p *Point3f) *Ray {
	origin := OffsetRayOrigin(i.point, i.pError, i.normal, p.Sub(i.point))
	direction := p.Sub(i.point)
	return &Ray{origin, direction, 1 - ShadowEpsilon, i.time, i.GetMedium(direction)}
}

func (i *Interaction) SpawnRayToInteraction(to *Interaction) *Ray {
	origin := OffsetRayOrigin(i.point, i.pError, i.normal, to.point.Sub(i.point))
	target := OffsetRayOrigin(to.point, to.pError, to.normal, origin.Sub(to.point))
	direction := target.Sub(origin)

	return &Ray{
		origin: origin,
		direction: direction,
		tMax: 1 - ShadowEpsilon,
		time: i.time,
		medium: i.GetMedium(direction),
	}
}

func (i *Interaction) IsSurfaceInteraction() bool {
	return i.normal != nil
}

func (i *Interaction) IsMediumInteraction() bool {
	return !i.IsSurfaceInteraction()
}

func (i *Interaction) GetMedium(w *Vector3f) Mediumer {
	if w.Dot(i.normal) > 0 {
		return i.mediumAccessor.Outside
	}
	return i.mediumAccessor.Inside
}

type PhaseFunction interface {
	P(wo, wi *Vector3f) float64
	SampleP(wo, wi *Vector3f, u *Point2f) float64
}

func PhaseHG(cosTheta float64, g float64) float64 {
	denom := 1.0 + g * g + 2 * g * cosTheta
	return Inv4Pi * (1.0 - g * g) / (denom * math.Sqrt(denom))
}

type Shading struct {
	normal     *Normal3f
	dpdu, dpdv *Vector3f
	dndu, dndv *Normal3f
}

type SurfaceInteraction struct {
	*Interaction

	uv                     *Point2f
	dpdu, dpdv             *Vector3f
	dndu, dndv             *Normal3f
	shape                  Shaper
	shading                *Shading
	primitive              Primitive
	bsdf                   *BSDF
	bssrdf                 *BSSRDF
	dpdx, dpdy             *Vector3f
	dudx, dvdx, dudy, dvdy float64

	// Shapes can optionally provide a face index with an
	// intersection point for use in Ptex texture lookups.
	// If Ptex isn't being used, then this value is ignored.
	faceIndex              int
}

func NewSurfaceInteraction(p *Point3f, pError *Vector3f, uv *Point2f, wo *Vector3f, dpdu, dpdv *Vector3f, dndu, dndv *Normal3f, time float64, shape Shaper, faceIndex int) *SurfaceInteraction {
	normal := Normal3f(*dpdu.Cross(dpdv).Normalized())

	return &SurfaceInteraction{
		Interaction: &Interaction{
			point: p,
			pError: pError,
			time: time,
			wo: wo,
			normal: &normal,
		},
		uv: uv,
		dpdu: dpdu,
		dpdv: dpdv,
		dndu: dndu,
		dndv: dndv,
		shape: shape,
		shading: &Shading{
			normal: &normal,
			dpdu:   dpdu,
			dpdv:   dpdv,
			dndu:   dndu,
			dndv:   dndv,
		},
		faceIndex: faceIndex,
	}
}

func (si *SurfaceInteraction) Le(w *Vector3f) Spectrum {
	area := si.primitive.GetAreaLight()
	if area != nil {
		return area.L(si, w)
	}
	return NewSpectrum(0.0)
}

type MediumInteraction struct {
	*Interaction

	phase *PhaseFunction
}

func (mi *MediumInteraction) IsValid() bool {
	return mi.phase != nil
}

type HenyeyGreenstein struct {
	g float64
}

func (hg *HenyeyGreenstein) P(wo, wi *Vector3f) float64 {
	return PhaseHG(wo.Dot(wi), hg.g)
}

func (hg *HenyeyGreenstein) SampleP(wo, wi *Vector3f, u *Point2f) float64 {
	var cosTheta float64
	if math.Abs(hg.g) < 1e-3 {
		cosTheta = 1.0 - 2.0 * u.X
	} else {
		sqrTerm := (1.0 - hg.g * hg.g) / (1.0 - hg.g + 2 * hg.g * u.X)
		cosTheta = (1.0 + hg.g * hg.g - sqrTerm * sqrTerm) / (2.0 * hg.g)
	}

	sinTheta := math.Sqrt(math.Max(0.0, 1.0 - cosTheta * cosTheta))
	phi := 2.0 * math.Pi * u.Y
	v1, v2 := CoordinateSystem(wo)
	*wi = *SphericalDirectionXYZ(sinTheta, cosTheta, phi, v1, v2, wo.MulScalar(-1.0))
	return PhaseHG(-cosTheta, hg.g)
}