package pbrt

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
)

type Distribution interface {
	D(wh *Vector3f) float64
	Lambda(w *Vector3f) float64
	G1(w *Vector3f) float64
	G(wo, wi *Vector3f) float64
	SampleWH(wo *Vector3f, u *Point2f) *Vector3f
	Pdf(wo, wh *Vector3f) float64

	SampleVisibleArea() bool
}

func G1(md Distribution, w *Vector3f) float64 {
	return 1 / (1 + md.Lambda(w))
}

func G(md Distribution, wo, wi *Vector3f) float64 {
	return 1 / (1 + md.Lambda(wo) + md.Lambda(wi))
}

func DistributionPdf(md Distribution, wo, wh *Vector3f) float64 {
	if md.SampleVisibleArea() {
		return md.D(wh) * md.G1(wo) * wo.AbsDot(wh) / AbsCosTheta(wo)
	} else {
		return md.D(wh) * AbsCosTheta(wh)
	}
}

func NewTrowbridgeReitz(alphaX, alphaY float64) Distribution {
	return &TrowbridgeReitz{
		alphaX:            alphaX,
		alphaY:            alphaY,
		sampleVisibleArea: true,
	}
}

type TrowbridgeReitz struct {
	alphaX, alphaY    float64
	sampleVisibleArea bool
}

func (t *TrowbridgeReitz) D(wh *Vector3f) float64 {
	tan2Theta := Tan2Theta(wh)
	if math.IsInf(tan2Theta) {
		return 0
	}
	cos4Theta := Cos2Theta(wh) * Cos2Theta(wh)
	e := (Cos2Phi(wh)/(t.alphaX*t.alphaY) + Sin2Phi(wh)/(t.alphaY*t.alphaY)) * tan2Theta
	return 1 / (math.Pi * t.alphaX * t.alphaY * cos4Theta * (1 + e) * (1 + e))
}
func (t *TrowbridgeReitz) Lambda(w *Vector3f) float64 {
	absTanTheta := math.Abs(TanTheta(w))
	if math.IsInf(absTanTheta) {
		return 0
	}

	alpha := math.Sqrt(Cos2Phi(w)*t.alphaX*t.alphaX + Sin2Phi(w)*t.alphaY*t.alphaY)
	alpha2Tan2Theta := (alpha * absTanTheta) * (alpha * absTanTheta)
	return (-1 + math.Sqrt(1.0+alpha2Tan2Theta)) / 2
}
func (t *TrowbridgeReitz) G1(w *Vector3f) float64 {
	return G1(t, w)
}
func (t *TrowbridgeReitz) G(wo, wi *Vector3f) float64 {
	return G(t, wo, wi)
}
func (t *TrowbridgeReitz) SampleWH(wo *Vector3f, u *Point2f) *Vector3f {
	var wh *Vector3f
	if !t.sampleVisibleArea {
		var cosTheta float64
		phi := math.Pi2 * u.Y
		if t.alphaX == t.alphaY {
			tanTheta2 := t.alphaX * t.alphaX * u.X / (1.0 - u.X)
			cosTheta = 1 / math.Sqrt(1+tanTheta2)
		} else {
			phi := math.Atan(t.alphaY / t.alphaX * math.Tan(math.Pi2*u.Y+0.5*math.Pi))
			if u.Y > 0.5 {
				phi += math.Pi
			}
			sinPhi := math.Sin(phi)
			cosPhi := math.Cos(phi)
			alphaX2 := t.alphaX * t.alphaX
			alphaY2 := t.alphaY * t.alphaY
			alpha2 := 1 / (cosPhi*cosPhi/alphaX2 + sinPhi*sinPhi/alphaY2)
			tanTheta2 := alpha2 * u.X / (1 - u.X)
			cosTheta = 1.0 / math.Sqrt(1+tanTheta2)
		}
		sinTheta := math.Sqrt(math.Max(0, 1-cosTheta*cosTheta))
		wh := SphericalDirection(sinTheta, cosTheta, phi)
		if !SameHemisphere(wo, wh) {
			wh = wh.MulScalar(-1)
		}
	} else {
		flip := wo.Z < 0
		if flip {
			wo = wo.MulScalar(-1)
		}
		wh := TrowbridgeReitzSample(wo, t.alphaX, t.alphaY, u.X, u.Y)
		if flip {
			wh = wh.MulScalar(-1)
		}
	}
	return wh
}
func (t *TrowbridgeReitz) Pdf(wo, wh *Vector3f) float64 {
	return DistributionPdf(t, wo, wh)
}

func (t *TrowbridgeReitz) SampleVisibleArea() bool {
	return t.sampleVisibleArea
}

func TrowbridgeReitzSample11(cosTheta, U1, U2 float64) (slopeX, slopeY float64) {
	// special case (normal incidence)
	if cosTheta > 0.9999 {
		r := math.Sqrt(U1 / (1 - U1))
		phi := math.Pi2 * U2
		return r * math.Cos(phi), r * math.Sin(phi)
	}

	sinTheta := math.Sqrt(math.Max(0, 1-cosTheta*cosTheta))
	tanTheta := sinTheta / cosTheta
	a := 1 / tanTheta
	G1 := 2 / (1 + math.Sqrt(1+1/(a*a)))

	// sample slopeX
	A := 2*U1/G1 - 1
	tmp := 1 / (A*A - 1)
	// TODO: is this necessary?
	if tmp > 1e10 {
		tmp = 1e10
	}

	B := tanTheta
	D := math.Sqrt(math.Max(B*B*tmp*tmp-(A*A-B*B)*tmp, 0))
	slopeX1 := B*tmp - D
	slopeX2 := B*tmp + D
	if A < 0 || slopeX2 > 1/tanTheta {
		slopeX = slopeX1
	} else {
		slopeX = slopeX2
	}

	// sample slopeY
	var S float64
	if U2 > 0.5 {
		S = 1.0
		U2 = 2.0 * (U2 - 0.5)
	} else {
		S = -1
		U2 = 2.0 * (0.5 - U2)
	}
	z := (U2 * (U2*(U2*0.27385-0.73369) + 0.46341)) /
		(U2*(U2*(U2*0.093073+0.309420)-1.000000) + 0.597999)
	slopeY = S * z * math.Sqrt(1+slopeX*slopeX)

	return slopeX, slopeY
}

func TrowbridgeReitzSample(wi *Vector3f, alphaX, alphaY, U1, U2 float64) *Vector3f {
	// 1. stretch wi
	wiStretched := &Vector3f{alphaX * wi.X, alphaY * wi.Y, wi.Z}
	wiStretched.Normalize()

	// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
	slopeX, slopeY := TrowbridgeReitzSample11(CosTheta(wiStretched), U1, U2)

	// 3. rotate
	slopeX, slopeY = CosPhi(wiStretched)*slopeX-SinPhi(wiStretched)*slopeY,
		SinPhi(wiStretched)*slopeX+CosPhi(wiStretched)*slopeY

	// 4. unstretch
	slopeX, slopeY = alphaX*slopeX, alphaY*slopeY

	// 5. compute normal
	v := &Vector3f{-slopeX, -slopeY, 1.}
	return v.Normalized()
}

// RoughnessToAlpha TrowbridgeReitzDistribution
func RoughnessToAlpha(roughness float64) float64 {
	roughness = math.Max(roughness, 1e-3)
	x := math.Log(roughness)
	return 1.62142 + 0.819955*x + 0.1734*x*x + 0.0171201*x*x*x + 0.000640711*x*x*x*x
}
