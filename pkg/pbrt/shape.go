package pbrt

import (
	"math"
)

type Shaper interface {
	ObjectBound() *Bounds3
	WorldBound() *Bounds3
	Intersect(r *Ray, testAlphaTexture bool) (intersects bool, tHit float64, si *SurfaceInteraction)
	IntersectP(r *Ray, testAlphaTexture bool) (intersects bool)
	Sample(u *Point2f) (i Interaction, pdf float64)
	SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64)
	Pdf(ref Interaction) float64
	PdfWi(ref Interaction) (float64, *Vector3f)
	SolidAngle(p *Point3f, nSamples int) float64
	Area() float64
}

type Shape struct {
	objectToWorld, worldToObject *Transform
	reverseOrientation           bool
	transformSwapsHandedness     bool
}

func (s *Shape) ObjectBound() *Bounds3 {
	return nil
}

func (s *Shape) WorldBound() *Bounds3 {
	return s.objectToWorld.TransformBounds(s.ObjectBound())
}

func (s *Shape) Intersect(r *Ray, testAlphaTexture bool) (intersects bool, tHit float64, si *SurfaceInteraction) {
	// TODO
	return false, 0, nil
}

func (s *Shape) IntersectP(r *Ray, testAlphaTexture bool) (intersects bool, tHit float64, si *SurfaceInteraction) {
	return s.Intersect(r, testAlphaTexture)
}

func (s *Shape) Area() float64 {
	// TODO
	return 0
}

func (s *Shape) Sample(u *Point2f) (i Interaction, pdf float64) {
	// TODO
	return nil, 0
}

// Pdf calculates the Probability Distribution Function
func (s *Shape) Pdf(ref Interaction) float64 {
	return 1 / s.Area()
}

func (s *Shape) PdfWi(ref Interaction) (float64, *Vector3f) {
	return s.Pdf(ref), new(Vector3f)
}

func (s *Shape) SampleAtInteraction(ref Interaction, u *Point2f) (i Interaction, pdf float64) {
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

func (s *Shape) PdfAtInteraction(ref Interaction, wi *Vector3f) float64 {
	r := ref.SpawnRay(wi)

	intersects, _, lightInteraction := s.Intersect(r, false)
	if !intersects {
		return 0.0
	}

	pdf := ref.GetPoint().Sub(lightInteraction.GetPoint()).LengthSquared() / lightInteraction.normal.AbsDot(wi.MulScalar(-1)) * s.Area()
	if math.IsInf(pdf, 1) {
		return 0.0
	}
	return pdf
}

func (s *Shape) SolidAngle(p *Point3f, nSamples int) float64 {
	ref := &interaction{
		point: p,
		wo: &Vector3f{0, 0, 1},
	}

	var solidAngle float64
	var i uint64
	for i = 1; i < uint64(nSamples + 1); i++ {
		u := &Point2f{RadicalInverse(0, i), RadicalInverse(1, i)}
		shapeInteraction, pdf := s.SampleAtInteraction(ref, u)
		if pdf <= 0 {
			continue
		}
		intersects, _, _ := s.IntersectP(&Ray{p, shapeInteraction.GetPoint().Sub(p), 0.999999, 0.0, nil}, true)
		if intersects {
			solidAngle += 1 / pdf
		}
	}
	return solidAngle / float64(nSamples)
}
