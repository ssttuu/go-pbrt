package lights

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type Distant struct {
	flags          pbrt.LightFlag
	samples        int32
	MediumAccessor *pbrt.MediumAccessor

	L           pbrt.Spectrum
	wLight      *pbrt.Vector3f
	worldCenter *pbrt.Point3f
	worldRadius float64
}

func NewDistant(lightToWorld *pbrt.Transform, L pbrt.Spectrum, wLight *pbrt.Vector3f) pbrt.Light {
	return &Distant{
		flags:   pbrt.LightFlagDeltaDirection,
		samples: 4,
		L:       L,
		wLight:  lightToWorld.TransformVector(wLight).Normalized(),
	}
}

func (d *Distant) GetFlags() pbrt.LightFlag {
	return d.flags
}

func (d *Distant) GetSamples() int32 {
	return d.samples
}

func (d *Distant) Preprocess(s pbrt.Scene) {
	d.worldCenter, d.worldRadius = s.WorldBound().BoundingSphere()
}

func (d *Distant) SampleLi(ref pbrt.Interaction, u *pbrt.Point2f) (s pbrt.Spectrum, wi *pbrt.Vector3f, pdf float64, vis *pbrt.VisibilityTester) {
	pOutside := d.wLight.MulScalar(2 * d.worldRadius)
	i := pbrt.NewInteraction(pOutside, new(pbrt.Vector3f), new(pbrt.Vector3f), new(pbrt.Vector3f), ref.GetTime(), d.MediumAccessor)
	return d.L, d.wLight, 1, pbrt.NewVisibilityTester(ref, i)
}

func (d *Distant) Power() pbrt.Spectrum {
	return d.L.MulScalar(math.Pi * d.worldRadius * d.worldRadius)
}

func (d *Distant) Le(r *pbrt.Ray) pbrt.Spectrum {
	return pbrt.Le(r)
}

func (d *Distant) PdfLi(ref pbrt.Interaction, wi *pbrt.Vector3f) float64 {
	return 0
}

func (d *Distant) SampleLe(u1, u2 *pbrt.Point2f, time float64) (s pbrt.Spectrum, r *pbrt.Ray, nLight *pbrt.Normal3f, pdfPos, pdfDir float64) {
	// choose point on disk oriented toward infinite light direction
	v1, v2 := pbrt.CoordinateSystem(d.wLight)
	cd := pbrt.ConcentricSampleDisk(u1)
	pDisk := d.worldCenter.Add(v1.MulScalar(cd.X).Add(v2.MulScalar(cd.Y)).MulScalar(d.worldRadius))

	// set ray origin and direction for infinite light ray
	r = pbrt.NewRay(pDisk.Add(d.wLight.MulScalar(d.worldRadius)), d.wLight.MulScalar(-1), time)

	return d.L, r, r.Direction, 1 / (math.Pi * d.worldRadius * d.worldRadius), 1
}

func (d *Distant) PdfLe(r *pbrt.Ray, nLight *pbrt.Normal3f) (pdfPos, pdfDir float64) {
	return 1 / (math.Pi * d.worldRadius * d.worldRadius), 0
}
