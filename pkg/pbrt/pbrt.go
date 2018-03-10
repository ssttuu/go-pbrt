package pbrt

import "math"

var Infinity = math.Inf(1)

const (
	ShadowEpsilon = 0.0001
	Pi            = 3.14159265358979323846
	InvPi         = 0.31830988618379067154
	Inv2Pi        = 0.15915494309189533577
	Inv4Pi        = 0.07957747154594766788
	PiOver2       = 1.57079632679489661923
	PiOver4       = 0.78539816339744830961
	Sqrt2         = 1.41421356237309504880

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

func FindInterval(size int, pred func(int) bool) int {
	first := 0
	l := size

	for l > 0 {
		half := l >> 1
		middle := first + half
		// bisect range based on value of pred at middle
		if pred(middle) {
			first = middle + 1
			l -= half + 1
		} else {
			l = half
		}
	}
	return int(Clamp(float64(first - 1), 0, float64(size - 2)))
}
