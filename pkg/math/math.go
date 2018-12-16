package math

import "math"

var Infinity = math.Inf(1)

var (
	Pi      = math.Pi
	Pi2     = Pi * 2
	InvPi   = 1.0 / Pi
	Inv2Pi  = 1.0 / (Pi * 2)
	Inv4Pi  = 1.0 / (Pi * 4)
	PiOver2 = Pi / 2.0
	PiOver4 = Pi / 4.0
	Sqrt2   = math.Sqrt2

	MachineEpsilon  = NextFloatUp(0.0)
	OneMinusEpsilon = NextFloatDown(1.0)
	ShadowEpsilon   = 0.0001
)

func Abs(x float64) float64 {
	return math.Abs(x)
}

func Acos(x float64) float64 {
	return math.Acos(x)
}

func Atan(x float64) float64 {
	return math.Atan(x)
}

func Atan2(y, x float64) float64 {
	return math.Atan2(y, x)
}

func Ceil(x float64) float64 {
	return math.Ceil(x)
}

func Clamp(v, low, high float64) float64 {
	if v < low {
		return low
	}
	if v > high {
		return high
	}
	return v
}

func Cos(x float64) float64 {
	return math.Cos(x)
}

func Degrees(rad float64) float64 {
	return 180.0 / Pi * rad
}

func Floor(x float64) float64 {
	return math.Floor(x)
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
	return int(Clamp(float64(first-1), 0, float64(size-2)))
}

func Gamma(n float64) float64 {
	return (n * MachineEpsilon) / (1 - n*MachineEpsilon)
}

func Inf(sign int) float64 {
	return math.Inf(sign)
}

func IsNegativeInf(f float64) bool {
	return math.IsInf(f, -1)
}

func IsPositiveInf(f float64) bool {
	return math.IsInf(f, 1)
}

func IsInf(f float64) bool {
	return math.IsInf(f, 1) || math.IsInf(f, -1)
}

func IsNaN(f float64) bool {
	return math.IsNaN(f)
}

func Lerp(t, v1, v2 float64) float64 {
	return (1.0-t)*v1 + t*v2
}

func Log(x float64) float64 {
	return math.Log(x)
}

func Max(x, y float64) float64 {
	return math.Max(x, y)
}

func Min(x, y float64) float64 {
	return math.Min(x, y)
}

func NextFloatUp(v float64) float64 {
	return math.Nextafter(v, v+1)
}

func NextFloatDown(v float64) float64 {
	return math.Nextafter(v, v-1)
}

func Radians(deg float64) float64 {
	return Pi / 180.0 * deg
}

func Sin(x float64) float64 {
	return math.Sin(x)
}

func Sqrt(x float64) float64 {
	return math.Sqrt(x)
}

func Tan(x float64) float64 {
	return math.Tan(x)
}
