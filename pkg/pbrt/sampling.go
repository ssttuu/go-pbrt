package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

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
		d.Cdf[i] = d.Cdf[i-1] + d.Func[i-1]/float64(n)
	}

	// transform step function integral into CDF
	d.FuncInt = d.Cdf[n]
	if d.FuncInt == 0.0 {
		for i := 1; i < n+1; i++ {
			d.Cdf[i] = float64(i) / float64(n)
		}
	} else {
		for i := 1; i < n+1; i++ {
			d.Cdf[i] /= d.FuncInt
		}
	}

	return d
}

func (d *Distribution1D) Count() int {
	return len(d.Func)
}

func (d *Distribution1D) SampleDiscrete(u float64) (value int, pdf float64, uRemapped float64) {
	offset := math.FindInterval(len(d.Cdf), func(index int) bool {
		return d.Cdf[index] <= u
	})

	pdf = 0
	if d.FuncInt > 0 {
		pdf = d.Func[offset] / (d.FuncInt / float64(d.Count()))
	}

	uRemapped = (u - d.Cdf[offset]) / (d.Cdf[offset+1] - d.Cdf[offset])

	return offset, pdf, uRemapped
}

func LatinHypercube(samples []float64, nSamples, nDim int32, rng *RandomNumberGenerator) {
	// Generate LHS samples along diagonal
	invNSamples := 1.0 / float64(nSamples)
	for i := int32(0); i < nSamples; i++ {
		for j := int32(0); j < nDim; j++ {
			sj := (float64(i) + rng.UniformFloat()) * invNSamples
			samples[nDim * i + j] = math.Min(sj, math.OneMinusEpsilon)
		}
	}

	// Permute LHS samples in each dimension
	for i := int32(0); i < nDim; i++ {
		for j := int32(0); j < nSamples; j++ {
			other := j + int32(rng.UniformUInt32B(uint32(nSamples - j)))
			samples[nDim * j + i], samples[nDim * other + i] = samples[nDim * other + i], samples[nDim * j + i]
		}
	}
}

func LatinHypercube2D(samples []Point2f, rng *RandomNumberGenerator) {
	// Generate LHS samples along diagonal
	nSamples := len(samples)
	invNSamples := 1.0 / float64(nSamples)
	for i := 0; i < nSamples; i++ {
		// for X and Y
		for j := 0; j < 2; j++ {
			sj := (float64(i) + rng.UniformFloat()) * invNSamples
			samples[i].SetIndex(j, math.Min(sj, math.OneMinusEpsilon))
		}
	}

	// Permute LHS samples in each dimension
	for i := 0; i < 2; i++ {
		for j := 0; j < nSamples; j++ {
			other := j + int(rng.UniformUInt32B(uint32(nSamples - j)))

			// swap
			tmp := samples[j].GetIndex(i)
			samples[j].SetIndex(i, samples[other].GetIndex(i))
			samples[other].SetIndex(i, tmp)
		}
	}
}

func StratifiedSample1D(samp []float64, nSamples int32, rng *RandomNumberGenerator, jitter bool) {
	invNSamples := 1.0 / float64(nSamples)
	for i := int32(0); i < nSamples; i++ {
		delta := 0.5
		if jitter {
			delta = rng.UniformFloat()
		}
		samp[i] = math.Min((float64(i)+delta) * invNSamples, math.OneMinusEpsilon)
	}
}

func StratifiedSample2D(samp []Point2f, nx, ny int32, rng *RandomNumberGenerator, jitter bool) {
	var dx float64 = 1.0 / float64(nx)
	var dy float64 = 1.0 / float64(ny)
	for y := int32(0); y < ny; y++ {
		for x := int32(0); x < nx; x++ {
			jx, jy := 0.5, 0.5
			if jitter {
				jx = rng.UniformFloat()
				jy = rng.UniformFloat()
			}
			s := samp[y * nx + x]
			s.X = math.Min((float64(x) + jx) * dx, math.OneMinusEpsilon)
			s.X = math.Min((float64(y) + jy) * dy, math.OneMinusEpsilon)
		}
	}
}

func ShuffleSamples1D(samp []float64, count, nDimensions int32, rng *RandomNumberGenerator) {
	for i := int32(0); i < count; i++ {
		other := i + int32(rng.UniformUInt32B(uint32(count - i)))
		for j := int32(0); j < nDimensions; j++ {
			samp[nDimensions * i + j], samp[nDimensions * other + j] = samp[nDimensions * other + j], samp[nDimensions * i + j]
		}
	}
}

func ShuffleSamples2D(samp []Point2f, count, nDimensions int32, rng *RandomNumberGenerator) {
	for i := int32(0); i < count; i++ {
		other := i + int32(rng.UniformUInt32B(uint32(count - i)))
		for j := int32(0); j < nDimensions; j++ {
			samp[nDimensions * i + j], samp[nDimensions * other + j] = samp[nDimensions * other + j], samp[nDimensions * i + j]
		}
	}
}

func UniformSampleHemisphere(u *Point2f) *Vector3f {
	z := u.X
	r := math.Sqrt(math.Max(0, 1.0-z*z))
	phi := 2 * math.Pi * u.Y
	return &Vector3f{r * math.Cos(phi), r * math.Sin(phi), z}
}

func UniformHemispherePdf() float64 {
	return math.Inv2Pi
}

func UniformSampleSphere(u *Point2f) *Vector3f {
	z := 1.0 - 2.0*u.X
	r := math.Sqrt(math.Max(0, 1-z*z))
	phi := 2 * math.Pi * u.Y
	return &Vector3f{r * math.Cos(phi), r * math.Sin(phi), z}
}

func UniformSpherePdf() float64 {
	return math.Inv4Pi
}

func UniformConePdf(cosThetaMax float64) float64 {
	return 1.0 / (2.0 * math.Pi * (1.0 - cosThetaMax))
}

func ConcentricSampleDisk(u *Point2f) *Point2f {
	uOffset := u.MulScalar(2.0).Sub(&Vector2f{1, 1})

	// handle degeneracy at the Origin
	if uOffset.X == 0 && uOffset.Y == 0 {
		return new(Point2f)
	}

	// apply concentric mapping to Point
	var theta, r float64
	if math.Abs(uOffset.X) > math.Abs(uOffset.Y) {
		r = uOffset.X
		theta = math.PiOver4 * (uOffset.Y / uOffset.X)
	} else {
		r = uOffset.Y
		theta = math.PiOver2 - math.PiOver4*(uOffset.X/uOffset.Y)
	}
	p := &Point2f{math.Cos(theta), math.Sin(theta)}
	return p.MulScalar(r)
}

func CosineSampleHemisphere(u *Point2f) *Vector3f {
	d := ConcentricSampleDisk(u)
	z := math.Sqrt(math.Max(0.0, 1.0-d.X*d.X-d.Y*d.Y))
	return &Vector3f{d.X, d.Y, z}
}

func CosineHemispherePdf(cosTheta float64) float64 {
	return cosTheta * math.InvPi
}

func BalanceHeuristic(nf int, fPdf float64, ng int, gPdf float64) float64 {
	return (float64(nf) * fPdf) / (float64(nf)*fPdf + float64(ng)*gPdf)
}

func PowerHeuristic(nf int, fPdf float64, ng int, gPdf float64) float64 {
	f := float64(nf) * fPdf
	g := float64(ng) * gPdf
	return (f * f) / (f*f + g*g)
}
