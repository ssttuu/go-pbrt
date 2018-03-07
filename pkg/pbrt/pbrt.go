package pbrt

import "math"

var Infinity = math.Inf(1)

const (
	Inv4Pi        = 0.07957747154594766788
	ShadowEpsilon = 0.0001

	OneMinusEpsilon = 0.99999999999999989
)

func Clamp(v, low, high float64) float64 {
	if v < low {
		return low
	}
	if v > high {
		return high
	}
	return v
}

func Lerp(t, v1, v2 float64) float64 {
	return (1.0-t)*v1 + t*v2
}

func Radians(deg float64) float64 {
	return math.Pi / 180.0 * deg
}

func Degrees(rad float64) float64 {
	return 180.0 / math.Pi * rad
}

func NextFloatUp(v float64) float64 {
	if math.IsInf(v, 1) {
		return v
	}
	if v == -0.0 {
		v = 0.0
	}

	ui := uint64(v)
	if (v >= 0) {
		ui++
	} else {
		ui--
	}

	return float64(ui)
}

func NextFloatDown(v float64) float64 {
	if math.IsInf(v, -1) {
		return v
	}
	if v == 0.0 {
		v = -0.0
	}

	ui := uint64(v)
	if (v > 0) {
		ui--
	} else {
		ui++
	}

	return float64(ui)
}
