//go:generate mockgen -source=shape.go -destination=shape.mock.go -package=pbrt

package pbrt

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
)

type Shape interface {
	ObjectBound() *Bounds3
	WorldBound() *Bounds3
	ReverseOrientation() bool
	TransformSwapsHandedness() bool
	Intersect(r *Ray, si *SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64)
	IntersectP(r *Ray, testAlphaTexture bool) (intersects bool)
	Sample(u *Point2f) (i Interaction, pdf float64)
	SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64)
	Pdf(ref Interaction) float64
	PdfWi(ref Interaction, wi *Vector3f) float64
	SolidAngle(p *Point3f, nSamples int) float64
	Area() float64
}

// Pdf calculates the Probability Distribution Function
func Pdf(s Shape, ref Interaction) float64 {
	return 1 / s.Area()
}

func PdfWi(s Shape, ref Interaction, wi *Vector3f) float64 {
	ray := ref.SpawnRay(wi)
	// Ignore any alpha textures used for trimming the shape when performing
	// this intersection. Hack for the "San Miguel" scene, where this is used
	// to make an invisible area light.
	intersectLight := NewSurfaceInteraction()
	intersects, _ := s.Intersect(ray, intersectLight, false)
	if !intersects {
		return 0
	}

	// convert light sample weight to solid angle measure
	pdf := ref.GetPoint().DistanceSquared(intersectLight.Point) / intersectLight.Normal.AbsDot(wi.MulScalar(-1).MulScalar(s.Area()))
	if math.IsInf(pdf) {
		return 0
	}

	return pdf
}

func SampleAtInteraction(s Shape, ref Interaction, u *Point2f) (intr Interaction, pdf float64) {
	intr, pdf = s.Sample(u)
	wi := intr.GetPoint().Sub(ref.GetPoint())
	if wi.LengthSquared() == 0.0 {
		pdf = 0
	} else {
		wi.Normalize()
		// Convert from area measure, as returned by the Sample() call above, to solid angle measure.
		pdf *= ref.GetPoint().DistanceSquared(intr.GetPoint()) / intr.GetNormal().AbsDot(wi.MulScalar(-1))
		if math.IsInf(pdf) {
			pdf = 0
		}
	}

	return intr, pdf

}

func SolidAngle(s Shape, p *Point3f, nSamples int) float64 {
	ref := NewInteraction(p, new(Vector3f), new(Normal3f), &Vector3f{0, 0, 1}, 0, &MediumAccessor{})

	solidAngle := 0.0
	for i := uint64(0); i < uint64(nSamples); i++ {
		u := &Point2f{X: RadicalInverse(0, i), Y: RadicalInverse(1, i)}
		pShape, pdf := s.SampleAtInteraction(ref, u)
		if pdf > 0 && !s.IntersectP(NewRay(p, pShape.GetPoint().Sub(p), 0.999), false) {
			solidAngle += 1.0 / pdf
		}
	}
	return solidAngle / float64(nSamples)
}
