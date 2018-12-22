package lights

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type DiffuseAreaLight struct {
	*pbrt.AreaLight

	LEmit    pbrt.Spectrum
	shape    pbrt.Shape
	twoSided bool
	area     float64
}

func NewDiffuseAreaLight(lightToWorld *pbrt.Transform, mediumAccessor *pbrt.MediumAccessor, LEmit pbrt.Spectrum, nSamples int32, shape pbrt.Shape, twoSided bool) pbrt.AreaLighter {
	return &DiffuseAreaLight{
		AreaLight: pbrt.NewAreaLight(lightToWorld, mediumAccessor, nSamples),
		LEmit:     LEmit,
		shape:     shape,
		twoSided:  twoSided,
		area:      shape.Area(),
	}
}

func (l *DiffuseAreaLight) Power() pbrt.Spectrum {
	multiplier := 1.0
	if l.twoSided {
		multiplier = 2
	}
	return l.LEmit.MulScalar(multiplier * l.area * math.Pi)

}

func (l *DiffuseAreaLight) L(intr pbrt.Interaction, w *pbrt.Vector3f) pbrt.Spectrum {
	if l.twoSided || intr.GetNormal().Dot(w) > 0 {
		return l.LEmit
	}
	return pbrt.NewSpectrum(0)
}

func (l *DiffuseAreaLight) Le(r *pbrt.Ray) pbrt.Spectrum {
	return pbrt.Le(r)
}

func (l *DiffuseAreaLight) SampleLi(ref pbrt.Interaction, u *pbrt.Point2f) (s pbrt.Spectrum, wi *pbrt.Vector3f, pdf float64, vis *pbrt.VisibilityTester) {
	// si := l.AreaLight.LightToWorld.TransformSurfaceInteraction(ref.(*pbrt.SurfaceInteraction))
	pShape, pdf := l.shape.SampleAtInteraction(ref, u)
	pShape.SetMediumAccessor(l.MediumAccessor)
	if pdf == 0 || pShape.GetPoint().Sub(ref.GetPoint()).LengthSquared() == 0 {
		return pbrt.NewSpectrum(0), new(pbrt.Vector3f), 0, new(pbrt.VisibilityTester)
	}

	wi = pShape.GetPoint().Sub(ref.GetPoint())
	vis = pbrt.NewVisibilityTester(ref, pShape)

	return l.L(pShape, wi.MulScalar(-1)), wi, pdf, vis
}

func (l *DiffuseAreaLight) PdfLi(ref pbrt.Interaction, wi *pbrt.Vector3f) float64 {
	return l.shape.PdfWi(ref, wi)
}

func (l *DiffuseAreaLight) SampleLe(u1, u2 *pbrt.Point2f, time float64) (s pbrt.Spectrum, r *pbrt.Ray, nLight *pbrt.Normal3f, pdfPos, pdfDir float64) {
	// sample a Point on the are light's shape, pShape
	pShape, pdfPos := l.shape.Sample(u1)
	pShape.SetMediumAccessor(l.MediumAccessor)

	var w *pbrt.Vector3f
	if l.twoSided {
		u := u2
		// choose a side to sample and then remap u[0] to [0,1] before
		// applying cosine-weighted hemispher sampling for the chose side.
		if u.X < 0.5 {
			u.X = math.Min(u.X*2, math.OneMinusEpsilon)
			w = pbrt.CosineSampleHemisphere(u)
		} else {
			u.X = math.Min((u.X-0.5)*2, math.OneMinusEpsilon)
			w = pbrt.CosineSampleHemisphere(u)
			w.Z *= -1
		}
		pdfDir = 0.5 * pbrt.CosineHemispherePdf(math.Abs(w.Z))
	} else {
		w = pbrt.CosineSampleHemisphere(u2)
		pdfDir = pbrt.CosineHemispherePdf(w.Z)
	}

	v1, v2 := pbrt.CoordinateSystem(pShape.GetNormal())
	w = v1.MulScalar(w.X).Add(v2.MulScalar(w.Y)).Add(pShape.GetNormal().MulScalar(w.Z))
	return l.L(pShape, w), pShape.SpawnRay(w), pShape.GetNormal(), pdfPos, pdfDir
}

func (l *DiffuseAreaLight) PdfLe(r *pbrt.Ray, n *pbrt.Normal3f) (pdfPos, pdfDir float64) {
	it := pbrt.NewInteraction(r.Origin, &pbrt.Vector3f{}, n, n, r.Time, l.MediumAccessor)
	pdfPos = l.shape.Pdf(it)
	if l.twoSided {
		pdfPos = 0.5 * pbrt.CosineHemispherePdf(n.AbsDot(r.Direction))
	} else {
		pdfPos = pbrt.CosineHemispherePdf(n.AbsDot(r.Direction))
	}
	return pdfPos, pdfDir
}
