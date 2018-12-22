package lights

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type Point struct {
	flags          pbrt.LightFlag
	nSamples       int32
	MediumAccessor *pbrt.MediumAccessor
	lightToWorld   *pbrt.Transform
	worldToLight   *pbrt.Transform

	pLight *pbrt.Point3f
	I      pbrt.Spectrum
}

func NewPoint(lightToWorld *pbrt.Transform, mediumAccessor *pbrt.MediumAccessor, i pbrt.Spectrum) pbrt.Light {
	p, _ := lightToWorld.TransformPoint(new(pbrt.Point3f), new(pbrt.Vector3f))
	return &Point{
		flags:          pbrt.LightFlagDeltaPosition,
		nSamples:       8,
		MediumAccessor: mediumAccessor,
		lightToWorld:   lightToWorld,
		worldToLight:   lightToWorld.Inverse(),
		pLight:         p,
		I:              i,
	}
}

func (l *Point) GetFlags() pbrt.LightFlag {
	return l.flags
}

func (l *Point) GetSamples() int32 {
	return l.nSamples
}

func (l *Point) Preprocess(s pbrt.Scene) {

}

func (l *Point) SampleLi(ref pbrt.Interaction, u *pbrt.Point2f) (s pbrt.Spectrum, wi *pbrt.Vector3f, pdf float64, vis *pbrt.VisibilityTester) {
	wi = l.pLight.Sub(ref.GetPoint()).Normalized()
	pdf = 1.0
	vis = pbrt.NewVisibilityTester(ref, pbrt.NewInteraction(l.pLight, new(pbrt.Vector3f), new(pbrt.Vector3f), new(pbrt.Vector3f), ref.GetTime(), l.MediumAccessor))
	return l.I.DivScalar(l.pLight.DistanceSquared(ref.GetPoint())), wi, pdf, vis
}

func (l *Point) Power() pbrt.Spectrum {
	return l.I.MulScalar(4 * math.Pi)
}

func (l *Point) Le(r *pbrt.Ray) pbrt.Spectrum {
	return pbrt.Le(r)
}

func (l *Point) PdfLi(ref pbrt.Interaction, wi *pbrt.Vector3f) float64 {
	return 0
}

func (l *Point) SampleLe(u1, u2 *pbrt.Point2f, time float64) (s pbrt.Spectrum, r *pbrt.Ray, nLight *pbrt.Normal3f, pdfPos, pdfDir float64) {
	r = pbrt.NewRayWithMedium(l.pLight, pbrt.UniformSampleSphere(u1), time, l.MediumAccessor.Inside)
	return l.I, r, r.Direction, 1, pbrt.UniformSpherePdf()
}

func (l *Point) PdfLe(r *pbrt.Ray, nLight *pbrt.Normal3f) (pdfPos, pdfDir float64) {
	return 0, pbrt.UniformSpherePdf()
}
