package efloat

import "github.com/stupschwartz/go-pbrt/pkg/math"

func Abs(f *EFloat) *EFloat {
	if f.Low >= 0 {
		// The entire interval is greater than zero, so we're all set.
		return f
	} else if f.High <= 0 {
		// The entire interval is less than zero
		r := f.MulScalar(-1.0)
		r.Check()
		return r
	} else {
		// The interval straddle zero
		r := f
		r.Value = math.Abs(r.Value)
		r.Low = 0.0
		r.High = math.Max(-f.Low, f.High)
		r.Check()
		return r
	}
}

func Sqrt(f *EFloat) *EFloat {
	r := &EFloat{
		Value: math.Sqrt(f.Value),
		Low:   math.NextFloatDown(math.Sqrt(f.Low)),
		High:  math.NextFloatUp(math.Sqrt(f.High)),
	}
	r.Check()
	return r
}

func Quadratic(a, b, c *EFloat) (t0, t1 *EFloat, ok bool) {
	// Find quadratic discriminant
	discriminant := b.Value*b.Value - 4.*a.Value*c.Value
	if discriminant < 0 {
		return nil, nil, false
	}
	rootDiscriminant := math.Sqrt(discriminant)

	floatRootDiscriminant := New(rootDiscriminant, math.MachineEpsilon*rootDiscriminant)

	// Compute quadratic _t_ values
	q := &EFloat{}
	if b.Value < 0 {
		q = b.Sub(floatRootDiscriminant).MulScalar(-0.5)
	} else {
		q = b.Add(floatRootDiscriminant).MulScalar(-0.5)
	}

	t0 = q.Div(a)
	t1 = c.Div(q)
	if t0.Value > t1.Value {
		t0, t1 = t1, t0
	}
	return t0, t1, true
}
