package pbrt


type Light struct {
	flags int
	nSamples int
	MediumAccessor *MediumAccessor
	lightToWorld, worldToLight Transform

	//SampleLi(ref Interactioner, u *Point2f, wi *Vector3f, pdf float64, vis VisibilityTester) Spectrum
	//Power() Spectrum
	//Preprocess(*Scene)
	//Le(r *RayDifferential) Spectrum
	//PdfLi(ref Interactioner, wi *Vector3f) float64
	//SampleLe(u1, u2 *Point2f, time float64, ray *Ray, nLight *Normal3f, pdfPos, pdfDir float64) Spectrum
	//PdfLe(r *Ray, nLight *Normal3f, pdfPos, pdfDir float64)
}

type VisibilityTester struct {
	p0, p1 *Interaction
}

func (vt *VisibilityTester) Unoccluded(s *Scene) bool {
	return !s.IntersectP(vt.p0.SpawnRayToInteraction(vt.p1))
}

//func (vt *VisibilityTester) Tr(scene *Scene, sampler *Sampler) Spectrum {
//
//}

type AreaLighter interface {
	L(i Interactioner, w *Vector3f) Spectrum
}

type AreaLight struct {
	*Light
}
