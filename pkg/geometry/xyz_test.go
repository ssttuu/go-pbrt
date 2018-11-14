package geometry

import (
	"testing"
	"github.com/stretchr/testify/assert"
)


func TestXYZFloat64_Abs(t *testing.T) {
	x := &XYZFloat64{-1.0, -2.0, -3.0}
	expected := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, expected, x.Abs())
}

func TestXYZFloat64_AbsDot(t *testing.T) {
	x := &XYZFloat64{-1.0, -2.0, -3.0}
	assert.Equal(t, 14.0, x.AbsDot(&XYZFloat64{1.0, 2.0, 3.0}))
}

func TestXYZFloat64_Add(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, &XYZFloat64{2.0, 4.0, 6.0}, x.Add(&XYZFloat64{1.0, 2.0, 3.0}))
}

func TestXYZFloat64_AddAssign(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	x.AddAssign(&XYZFloat64{1.0, 2.0, 3.0})
	assert.Equal(t, &XYZFloat64{2.0, 4.0, 6.0}, x)
}

func TestXYZFloat64_AddConst(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, &XYZFloat64{2.0, 3.0, 4.0}, x.AddConst(1.0))
}

func TestXYZFloat64_Cross(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, &XYZFloat64{2.0, -1, 0.0}, x.Cross(&XYZFloat64{1.0, 2.0, 4.0}))
}

func TestXYZFloat64_Distance(t *testing.T) {
	x := &XYZFloat64{1, 2, 3}
	assert.Equal(t, 1.0, x.Distance(&XYZFloat64{1, 2, 4}))
}

func TestXYZFloat64_DistanceSquared(t *testing.T) {
	x := &XYZFloat64{1, 2, 3}
	assert.Equal(t, 4.0, x.DistanceSquared(&XYZFloat64{1, 2, 5}))
}

func TestXYZFloat64_Div(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	assert.Equal(t, &XYZFloat64{1, 3, 9}, x.Div(&XYZFloat64{3, 3, 3}))
}

func TestXYZFloat64_DivAssign(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	x.DivAssign(&XYZFloat64{3, 3, 3})
	assert.Equal(t, &XYZFloat64{1, 3, 9}, x)
}

func TestXYZFloat64_DivScalar(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	assert.Equal(t, &XYZFloat64{1, 3, 9}, x.DivScalar(3.0))
}

func TestXYZFloat64_Dot(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	y := &XYZFloat64{2, 4, 6}
	assert.Equal(t, 204.0, x.Dot(y))
}

func TestXYZFloat64_Equals(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	assert.True(t, x.Equals(x))
	assert.False(t, x.Equals(&XYZFloat64{}))
}

func TestXYZFloat64_Index(t *testing.T) {
	x := &XYZFloat64{2.0, 3.0, 4.0}
	assert.Equal(t, 2.0, x.Index(0))
	assert.Equal(t, 3.0, x.Index(1))
	assert.Equal(t, 4.0, x.Index(2))
}

func TestXYZFloat64_Length(t *testing.T) {
	x := &XYZFloat64{0.0, 3.0, 0.0}
	assert.Equal(t, 3.0, x.Length())
}

func TestXYZFloat64_LengthSquared(t *testing.T) {
	x := &XYZFloat64{0.0, 3.0, 0.0}
	assert.Equal(t, 9.0, x.LengthSquared())
}

func TestXYZFloat64_Mul(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	other := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, &XYZFloat64{1.0, 4.0, 9.0}, x.Mul(other))
}

func TestXYZFloat64_MulAssign(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	x.MulAssign(&XYZFloat64{1.0, 2.0, 3.0})
	assert.Equal(t, &XYZFloat64{1.0, 4.0, 9.0}, x)
}

func TestXYZFloat64_MulScalar(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, &XYZFloat64{2.0, 4.0, 6.0}, x.MulScalar(2.0))
}

func TestXYZFloat64_Normalize(t *testing.T) {
	x := &XYZFloat64{0.0, 0.0, 3.0}
	x.Normalize()
	assert.Equal(t, &XYZFloat64{0, 0, 1}, x)
}

func TestXYZFloat64_Normalized(t *testing.T) {
	x := &XYZFloat64{0.0, 0.0, 3.0}
	assert.Equal(t, &XYZFloat64{0, 0, 1}, x.Normalized())
}

func TestXYZFloat64_NotEquals(t *testing.T) {
	x := &XYZFloat64{3, 9, 27}
	assert.True(t, x.NotEquals(&XYZFloat64{}))
	assert.False(t, x.NotEquals(x))
}

func TestXYZFloat64_SetIndex(t *testing.T) {
	x := XYZFloat64{2.0, 3.0, 4.0}

	x.SetIndex(0, 2.1)
	x.SetIndex(1, 3.1)
	x.SetIndex(2, 4.1)

	assert.Equal(t, 2.1, x.Index(0))
	assert.Equal(t, 3.1, x.Index(1))
	assert.Equal(t, 4.1, x.Index(2))
}

func TestXYZFloat64_Set(t *testing.T) {
	x := XYZFloat64{2.0, 3.0, 4.0}
	x.Set(1.1, 2.1, 3.1)
	assert.Equal(t, XYZFloat64{1.1, 2.1, 3.1}, x)
}

func TestXYZFloat64_String(t *testing.T) {
	x := &XYZFloat64{1.0, 2.0, 3.0}
	assert.Equal(t, "XYZfloat64(1, 2, 3)", x.String())
}

func TestXYZFloat64_Sub(t *testing.T) {
	x := &XYZFloat64{2, 4, 6}
	assert.Equal(t, &XYZFloat64{1, 1, 1}, x.Sub(&XYZFloat64{1, 3, 5}))
}

func TestXYZFloat64_SubAssign(t *testing.T) {
	x := &XYZFloat64{2, 4, 6}
	x.SubAssign(&XYZFloat64{1, 3, 5})
	assert.Equal(t, &XYZFloat64{1, 1, 1}, x)
}
