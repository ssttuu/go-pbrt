package pbrt

import "math"

type Distribution1D struct {
	Func, Cdf []float64
	FuncInt   float64
}

func NewDistribution1D(f []float64) *Distribution1D {
	n := len(f)
	d := &Distribution1D{
		Func: f,
		Cdf:  make([]float64, n+1),
	}

	// compute integral of step function at xi
	d.Cdf[0] = 0
	for i := 1; i < n+1; i++ {
		d.Cdf[i] = d.Cdf[i - 1] + d.Func[i - 1] / float64(n)
	}

	// transform step function integral into CDF
	d.FuncInt = d.Cdf[n]
	if d.FuncInt == 0.0 {
		for i := 1; i < n + 1; i++ {
			d.Cdf[i] = float64(i) / float64(n)
		}
	} else {
		for i := 1; i < n + 1; i++ {
			d.Cdf[i] /= d.FuncInt
		}
	}

	return d
}

func (d *Distribution1D) Count() int {
	return len(d.Func)
}

func (d *Distribution1D) SampleDiscrete(u float64) (value int, pdf float64, uRemapped float64) {
	offset := FindInterval(len(d.Cdf), func (index int) bool {
		return d.Cdf[index] <= u
	})

	pdf = 0
	if d.FuncInt > 0 {
		pdf = d.Func[offset] / (d.FuncInt / float64(d.Count()))
	}

	uRemapped = (u - d.Cdf[offset]) / (d.Cdf[offset + 1] - d.Cdf[offset])

	return offset, pdf, uRemapped

}

func UniformSampleHemisphere(u *Point2f) *Vector3f {
	z := u.X
	r := math.Sqrt(math.Max(0, 1.0 - z * z))
	phi := 2 * math.Pi * u.Y
	return &Vector3f{r * math.Cos(phi), r * math.Sin(phi), z}
}

func UniformHemispherePdf() float64 {
	return Inv2Pi
}

func UniformSampleSphere(u *Point2f) *Vector3f {
	z := 1.0 - 2.0*u.X
	r := math.Sqrt(math.Max(0, 1-z*z))
	phi := 2 * math.Pi * u.Y
	return &Vector3f{r * math.Cos(phi), r * math.Sin(phi), z}
}

func UniformSpherePdf() float64 {
	return Inv4Pi
}

func UniformConePdf(cosThetaMax float64) float64 {
	return 1.0 / (2.0 * math.Pi * (1.0 - cosThetaMax))
}

func ConcentricSampleDisk(u *Point2f) *Point2f {
	uOffset := u.MulScalar(2.0).Sub(&Vector2f{1, 1})

	// handle degeneracy at the Origin
	if uOffset.X == 0 && uOffset.Y == 0 {
		return &Point2f{}
	}

	// apply concentric mapping to Point
	var theta, r float64
	if math.Abs(uOffset.X) > math.Abs(uOffset.Y) {
		r = uOffset.X
		theta = PiOver4 * (uOffset.Y / uOffset.X)
	} else {
		r = uOffset.Y
		theta = PiOver2 - PiOver4*(uOffset.X/uOffset.Y)
	}
	return new(Point2f).Set(math.Cos(theta), math.Sin(theta)).MulScalar(r)
}

func CosineSampleHemisphere(u *Point2f) *Vector3f {
	d := ConcentricSampleDisk(u)
	z := math.Sqrt(math.Max(0.0, 1.0-d.LengthSquared()))
	return &Vector3f{d.X, d.Y, z}
}

func CosineHemispherePdf(cosTheta float64) float64 {
	return cosTheta * InvPi
}

func BalanceHeuristic(nf int, fPdf float64, ng int, gPdf float64) float64 {
	return (float64(nf) * fPdf) / (float64(nf)*fPdf + float64(ng)*gPdf)
}

func PowerHeuristic(nf int, fPdf float64, ng int, gPdf float64) float64 {
	f := float64(nf) * fPdf
	g := float64(ng) * gPdf
	return (f * f) / (f*f + g*g)
}
