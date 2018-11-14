//go:generate mockgen -source=shape.go -destination=shape.mock.go -package=pbrt

package pbrt

import (
	"math"
)

type Shape interface {
	GetName() string
	ObjectBound() *Bounds3
	WorldBound() *Bounds3
	Intersect(r *Ray, si *SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64)
	IntersectP(r *Ray, testAlphaTexture bool) (intersects bool)
	Sample(u *Point2f) (i Interaction, pdf float64)
	SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64)
	Pdf(ref Interaction) float64
	PdfWi(ref Interaction) (float64, *Vector3f)
	SolidAngle(p *Point3f, nSamples int) float64
	Area() float64
}

type shape struct {
	name string

	objectToWorld, worldToObject *Transform
	reverseOrientation           bool
	transformSwapsHandedness     bool
}

func (s *shape) GetName() string {
	return s.name
}

func (s *shape) ObjectBound() *Bounds3 {
	return nil
}

func (s *shape) WorldBound() *Bounds3 {
	return s.objectToWorld.TransformBounds(s.ObjectBound())
}

func (s *shape) Intersect(r *Ray, si *SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64) {
	// TODO
	return false, 0
}

func (s *shape) IntersectP(r *Ray, si *SurfaceInteraction, testAlphaTexture bool) (intersects bool, tHit float64) {
	return s.Intersect(r, si, testAlphaTexture)
}

func (s *shape) Area() float64 {
	// TODO
	return 0
}

func (s *shape) Sample(u *Point2f) (i Interaction, pdf float64) {
	// TODO
	return nil, 0
}

// Pdf calculates the Probability Distribution Function
func (s *shape) Pdf(ref Interaction) float64 {
	return 1 / s.Area()
}

func (s *shape) PdfWi(ref Interaction) (float64, *Vector3f) {
	return s.Pdf(ref), new(Vector3f)
}

func (s *shape) SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64) {
	intr, pdf := s.Sample(u)
	wi := intr.GetPoint().Sub(ref.GetPoint())
	if wi.LengthSquared() == 0.0 {
		return intr, 0
	}
	wi.Normalize()
	pdf = ref.GetPoint().Sub(intr.GetPoint()).LengthSquared() / intr.GetNormal().AbsDot(wi.MulScalar(-1))
	if math.IsInf(pdf, 1) {
		pdf = 0
	}

	return intr, pdf

}

func (s *shape) PdfAtInteraction(ref Interaction, wi *Vector3f) float64 {
	r := ref.SpawnRay(wi)

	lightInteraction := NewSurfaceInteraction()
	intersects, _ := s.Intersect(r, lightInteraction, false)
	if !intersects {
		return 0.0
	}

	pdf := ref.GetPoint().Sub(lightInteraction.GetPoint()).LengthSquared() / lightInteraction.Normal.AbsDot(wi.MulScalar(-1)) * s.Area()
	if math.IsInf(pdf, 1) {
		return 0.0
	}
	return pdf
}

func (s *shape) SolidAngle(p *Point3f, nSamples int) float64 {
	ref := &interaction{
		Point: p,
		wo:    &Vector3f{0, 0, 1},
	}

	var solidAngle float64
	var i uint64
	for i = 1; i < uint64(nSamples + 1); i++ {
		u := &Point2f{RadicalInverse(0, i), RadicalInverse(1, i)}
		shapeInteraction, pdf := s.SampleAtInteraction(ref, u)
		if pdf <= 0 {
			continue
		}
		ray := Ray{p, shapeInteraction.GetPoint().Sub(p), 0.999999, 0.0, nil}
		si := NewSurfaceInteraction()
		intersects, _ := s.IntersectP(&ray, si,true)
		if intersects {
			solidAngle += 1 / pdf
		}
	}
	return solidAngle / float64(nSamples)
}
