package pbrt

import "math"

func UniformSampleSphere(u *Point2f) *Vector3f {
	z := 1.0 - 2.0*u.X
	r := math.Sqrt(math.Max(0, 1-z*z))
	phi := 2 * math.Pi * u.Y
	return &Vector3f{r * math.Cos(phi), r * math.Sin(phi), z}
}

func UniformConePdf(cosThetaMax float64) float64 {
	return 1.0 / (2.0 * math.Pi * (1.0 - cosThetaMax))
}
