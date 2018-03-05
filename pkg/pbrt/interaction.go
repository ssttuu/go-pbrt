package pbrt

import (
	"math"
)

type Interaction interface {
	//IsSurfaceInteraction() bool
	//IsMediumInteraction() bool
	//
	//SpawnRay(direction *XYZGetter) *Ray
	GetPoint() *Point3
	GetPointError() *Vector3
	GetTime() float64
	GetNormal() *Normal3
	GetMedium(*Vector3) *Medium
}


func SpawnRayTo(from, to Interaction) *Ray {
	origin := OffsetRayOrigin(from.GetPoint(), from.GetPointError(), from.GetNormal(), to.GetPoint().Sub(from.GetPoint()) )
	target := OffsetRayOrigin(to.GetPoint(), to.GetPointError(), to.GetNormal(), origin.Sub(to.GetPoint()))
	direction := target.Sub(origin)

	return &Ray{
		origin: origin,
		direction: direction,
		tMax: 1 - ShadowEpsilon,
		time: from.GetTime(),
		medium: from.GetMedium(direction),
	}
}

type PhaseFunction interface {
	P(wo, wi *Vector3) float64
	SampleP(wo, wi *Vector3, u *Point2) float64
}

func PhaseHG(cosTheta float64, g float64) float64 {
	denom := 1.0 + g * g + 2 * g * cosTheta
	return Inv4Pi * (1.0 - g * g) / (denom * math.Sqrt(denom))
}

type SurfaceInteraction struct {
	point *Point3
	time float64
	pError *Vector3
	wo *Vector3
	normal *Normal3
	medium Medium
}

func (si *SurfaceInteraction) IsSurfaceInteraction() bool {
	return true
}

func (si *SurfaceInteraction) IsMediumInteraction() bool {
	return false
}

type MediumInteraction struct {
	point *Point3
	time float64
	pError *Vector3
	wo *Vector3
	normal *Normal3
	medium Medium

	phase *PhaseFunction
}

func (si *MediumInteraction) IsSurfaceInteraction() bool {
	return false
}

func (si *MediumInteraction) IsMediumInteraction() bool {
	return true
}


type HenyeyGreenstein struct {
	g float64
}

func (hg *HenyeyGreenstein) P(wo, wi *Vector3) float64 {
	return PhaseHG(wo.Dot(wi), hg.g)
}

func (hg *HenyeyGreenstein) SampleP(wo, wi *Vector3, u *Point2) float64 {
	var cosTheta float64
	if math.Abs(hg.g) < 1e-3 {
		cosTheta = 1.0 - 2.0 * u.x
	} else {
		sqrTerm := (1.0 - hg.g * hg.g) / (1.0 - hg.g + 2 * hg.g * u.x)
		cosTheta = (1.0 + hg.g * hg.g - sqrTerm * sqrTerm) / (2.0 * hg.g)
	}

	sinTheta := math.Sqrt(math.Max(0.0, 1.0 - cosTheta * cosTheta))
	phi := 2.0 * math.Pi * u.y
	var v1, v2 Vector3
	CoordinateSystem(wo, &v1, &v2)
	*wi = *SphericalDirectionXYZ(sinTheta, cosTheta, phi, &v1, &v2, wo.MulScalar(-1.0))
	return PhaseHG(-cosTheta, hg.g)
}