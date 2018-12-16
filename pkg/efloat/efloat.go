package efloat

import "C"
import (
	"log"

	"github.com/ssttuu/go-pbrt/pkg/math"
)

func New(v, err float64) *EFloat {
	f := &EFloat{
		Value: v,
		Low:   v,
		High:  v,
	}
	if err != 0 {
		f.Low = math.NextFloatDown(v - err)
		f.High = math.NextFloatUp(v + err)
	}
	f.Check()
	return f
}

type EFloat struct {
	Value float64
	Low   float64
	High  float64
}

func (f *EFloat) Add(other *EFloat) *EFloat {
	r := *f
	r.AddAssign(other)
	return &r
}
func (f *EFloat) AddAssign(other *EFloat) {
	f.Value = f.Value + other.Value
	f.Low = math.NextFloatDown(f.Low + other.Low)
	f.High = math.NextFloatUp(f.High + other.High)
	f.Check()
}

func (f *EFloat) Div(other *EFloat) *EFloat {
	r := *f
	r.DivAssign(other)
	return &r
}
func (f *EFloat) DivAssign(other *EFloat) {
	f.Value = f.Value / other.Value
	if other.Low < 0 && other.High > 0 {
		f.Low = -math.Infinity
		f.High = math.Infinity
	} else {
		div := [4]float64{
			f.Low / other.Low,
			f.High / other.Low,
			f.Low / other.High,
			f.High / other.High,
		}
		f.Low = math.NextFloatDown(math.Min(math.Min(div[0], div[1]), math.Min(div[2], div[3])))
		f.High = math.NextFloatUp(math.Max(math.Max(div[0], div[1]), math.Max(div[2], div[3])))
	}

	f.Check()
}

func (f *EFloat) Mul(other *EFloat) *EFloat {
	r := *f
	r.MulAssign(other)
	return &r
}
func (f *EFloat) MulAssign(other *EFloat) {
	prod := [4]float64{
		f.Low * other.Low,
		f.High * other.Low,
		f.Low * other.High,
		f.High * other.High,
	}

	f.Value = f.Value * other.Value
	f.Low = math.NextFloatDown(math.Min(math.Min(prod[0], prod[1]), math.Min(prod[2], prod[3])))
	f.High = math.NextFloatUp(math.Max(math.Max(prod[0], prod[1]), math.Max(prod[2], prod[3])))

	f.Check()
}

func (f *EFloat) MulScalar(other float64) *EFloat {
	return f.Mul(New(other, 0.0))
}

func (f *EFloat) Sub(other *EFloat) *EFloat {
	r := *f
	r.SubAssign(other)
	return &r
}
func (f *EFloat) SubAssign(other *EFloat) {
	f.Value = f.Value - other.Value
	f.Low = math.NextFloatDown(f.Low - other.High)
	f.High = math.NextFloatUp(f.High - other.Low)
	f.Check()
}

func (f *EFloat) Check() {
	// TODO: return error?
	if math.IsInf(f.Low) || math.IsNaN(f.Low) || math.IsInf(f.High) || math.IsNaN(f.High) {
		log.Panicf("EFloat IsInf or IsNan: %+v\n", f)
	}

	if f.Low > f.High {
		log.Panicf("EFloat Low is greater than High: %+v\n", f)
	}
}

func (f *EFloat) GetAbsoluteError() float64 {
	return f.High - f.Low
}
