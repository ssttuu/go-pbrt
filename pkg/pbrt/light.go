package pbrt


type Light interface {
	SampleLi(ref Interaction, u *Point2, wi *Vector3, pdf float64, vis VisibilityTester) Spectrum
	Power() Spectrum
	Preprocess(*Scene)
	Le(r *RayDifferential) Spectrum
	PdfLi(ref Interaction, wi *Vector3) float64
	SampleLe(u1, u2 *Point2, time float64, ray *Ray, nLight *Normal3, pdfPos, pdfDir float64) Spectrum
	PdfLe(r *Ray, nLight *Normal3, pdfPos, pdfDir float64)
}

type VisibilityTester struct {
	p0, p1 Interaction
}

func (vt *VisibilityTester) Unoccluded(s *Scene) bool {
	return !s.IntersectP(SpawnRayTo(vt.p0, vt.p1))
}

//func (vt *VisibilityTester) Tr(scene *Scene, sampler *Sampler) Spectrum {
//
//}

type AreaLight struct {
	flags int
	nSamples int
	medium Medium

	lightToWorld, worldToLight Transform
}
