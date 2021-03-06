//go:generate genny -in=$GOFILE -out=../../pkg/geometry/$GOFILE gen "GenericType=int64,float64"
package geometry

import (
	"fmt"
	"github.com/cheekybits/genny/generic"
	"math"
)

type GenericType generic.Type

type XYGenericType struct {
	X, Y GenericType
}

func (xy *XYGenericType) GetIndex(i int) GenericType {
	switch i {
	case 0:
		return xy.X
	default:
		return xy.Y
	}
}

func (xy *XYGenericType) SetIndex(i int, v GenericType) {
	switch i {
	case 0:
		xy.X = v
	default:
		xy.Y = v
	}
}

func (xy *XYGenericType) Set(x, y GenericType) *XYGenericType {
	xy.X = x
	xy.Y = y
	return xy
}

func (xy *XYGenericType) LengthSquared() GenericType {
	return xy.X*xy.X + xy.Y*xy.Y
}

func (xy *XYGenericType) Length() GenericType {
	return GenericType(math.Sqrt(float64(xy.LengthSquared())))
}

func (xy *XYGenericType) String() string {
	return fmt.Sprintf("XYGenericType(%v, %v)", xy.X, xy.Y)
}

func (xy *XYGenericType) Add(other *XYGenericType) *XYGenericType {
	return &XYGenericType{xy.X + other.X, xy.Y + other.Y}
}

func (xy *XYGenericType) AddAssign(other *XYGenericType) {
	xy.X += other.X
	xy.Y += other.Y
}

func (xy *XYGenericType) AddConst(other GenericType) *XYGenericType {
	return &XYGenericType{xy.X + other, xy.Y + other}
}

func (xy *XYGenericType) Sub(other *XYGenericType) *XYGenericType {
	return &XYGenericType{xy.X - other.X, xy.Y - other.Y}
}

func (xy *XYGenericType) SubAssign(other *XYGenericType) {
	xy.X -= other.X
	xy.Y -= other.Y
}

func (xy *XYGenericType) Mul(other *XYGenericType) *XYGenericType {
	return &XYGenericType{xy.X * other.X, xy.Y * other.Y}
}

func (xy *XYGenericType) MulAssign(other *XYGenericType) {
	xy.X *= other.X
	xy.Y *= other.Y
}

func (xy *XYGenericType) MulScalar(other GenericType) *XYGenericType {
	return &XYGenericType{xy.X * other, xy.Y * other}
}

func (xy *XYGenericType) Div(other *XYGenericType) *XYGenericType {
	return &XYGenericType{xy.X / other.X, xy.Y / other.Y}
}

func (xy *XYGenericType) DivAssign(other *XYGenericType) {
	xy.X /= other.X
	xy.Y /= other.Y
}

func (xy *XYGenericType) DivScalar(other GenericType) *XYGenericType {
	return &XYGenericType{xy.X / other, xy.Y / other}
}

func (xy *XYGenericType) Equals(other *XYGenericType) bool {
	return xy.X == other.X && xy.Y == other.Y
}

func (xy *XYGenericType) NotEquals(other *XYGenericType) bool {
	return xy.X != other.X || xy.Y != other.Y
}

func (xy *XYGenericType) Floor() *XYGenericType {
	return &XYGenericType{
		X: GenericType(math.Floor(float64(xy.X))),
		Y: GenericType(math.Floor(float64(xy.Y))),
	}
}

func (xy *XYGenericType) Ceil() *XYGenericType {
	return &XYGenericType{
		X: GenericType(math.Ceil(float64(xy.X))),
		Y: GenericType(math.Ceil(float64(xy.Y))),
	}
}

type XYZGenericType struct {
	X, Y, Z GenericType
}

func (xyz *XYZGenericType) Index(i int) GenericType {
	switch i {
	case 0:
		return xyz.X
	case 1:
		return xyz.Y
	default:
		return xyz.Z
	}
}

func (xyz *XYZGenericType) SetIndex(i int, v GenericType) {
	switch i {
	case 0:
		xyz.X = v
	case 1:
		xyz.Y = v
	default:
		xyz.Z = v
	}
}

func (xyz *XYZGenericType) Set(x, y, z GenericType) *XYZGenericType {
	xyz.X = x
	xyz.Y = y
	xyz.Z = z
	return xyz
}

func (xyz *XYZGenericType) LengthSquared() GenericType {
	return xyz.X*xyz.X + xyz.Y*xyz.Y + xyz.Z*xyz.Z
}

func (xyz *XYZGenericType) Length() GenericType {
	return GenericType(math.Sqrt(float64(xyz.LengthSquared())))
}

func (xyz *XYZGenericType) String() string {
	return fmt.Sprintf("XYZGenericType(%v, %v, %v)", xyz.X, xyz.Y, xyz.Z)
}

func (xyz *XYZGenericType) Add(other *XYZGenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X + other.X,
		xyz.Y + other.Y,
		xyz.Z + other.Z,
	}
}

func (xyz *XYZGenericType) AddAssign(other *XYZGenericType) {
	xyz.X += other.X
	xyz.Y += other.Y
	xyz.Z += other.Z
}

func (xyz *XYZGenericType) AddConst(other GenericType) *XYZGenericType {
	return &XYZGenericType{xyz.X + other, xyz.Y + other, xyz.Z + other}
}

func (xyz *XYZGenericType) Sub(other *XYZGenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X - other.X,
		xyz.Y - other.Y,
		xyz.Z - other.Z,
	}
}

func (xyz *XYZGenericType) SubAssign(other *XYZGenericType) {
	xyz.X -= other.X
	xyz.Y -= other.Y
	xyz.Z -= other.Z
}

func (xyz *XYZGenericType) Mul(other *XYZGenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X * other.X,
		xyz.Y * other.Y,
		xyz.Z * other.Z,
	}
}

func (xyz *XYZGenericType) MulAssign(other *XYZGenericType) {
	xyz.X *= other.X
	xyz.Y *= other.Y
	xyz.Z *= other.Z
}

func (xyz *XYZGenericType) MulScalar(other GenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X * other,
		xyz.Y * other,
		xyz.Z * other,
	}
}

func (xyz *XYZGenericType) Div(other *XYZGenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X / other.X,
		xyz.Y / other.Y,
		xyz.Z / other.Z,
	}
}

func (xyz *XYZGenericType) DivAssign(other *XYZGenericType) {
	xyz.X /= other.X
	xyz.Y /= other.Y
	xyz.Z /= other.Z
}

func (xyz *XYZGenericType) DivScalar(other GenericType) *XYZGenericType {
	return &XYZGenericType{
		xyz.X / other,
		xyz.Y / other,
		xyz.Z / other,
	}
}

func (xyz *XYZGenericType) IsNegative() [3]int {
	out := [3]int{}
	if xyz.X < 0 {
		out[0] = 1
	}
	if xyz.Y < 0 {
		out[1] = 1
	}
	if xyz.Z < 0 {
		out[2] = 1
	}
	return out
}

func (xyz *XYZGenericType) Abs() *XYZGenericType {
	return &XYZGenericType{
		GenericType(math.Abs(float64(xyz.X))),
		GenericType(math.Abs(float64(xyz.Y))),
		GenericType(math.Abs(float64(xyz.Z))),
	}
}

func (xyz *XYZGenericType) Dot(other *XYZGenericType) GenericType {
	return xyz.X*other.X + xyz.Y*other.Y + xyz.Z*other.Z
}

func (xyz *XYZGenericType) AbsDot(other *XYZGenericType) GenericType {
	return GenericType(math.Abs(float64(xyz.Dot(other))))
}

func (xyz *XYZGenericType) Distance(other *XYZGenericType) GenericType {
	return GenericType(math.Sqrt(xyz.DistanceSquared(other)))
}

func (xyz *XYZGenericType) DistanceSquared(other *XYZGenericType) float64 {
	return float64(other.Sub(xyz).LengthSquared())
}

func (xyz *XYZGenericType) Cross(other *XYZGenericType) *XYZGenericType {
	return &XYZGenericType{(xyz.Y * other.Z) - (xyz.Z * other.Y), (xyz.Z * other.X) - (xyz.X * other.Z), (xyz.X * other.Y) - (xyz.Y * other.X)}
}

func (xyz *XYZGenericType) Normalize() {
	nor2 := xyz.LengthSquared()
	if nor2 > 0 {
		invNor := GenericType(1) / GenericType(math.Sqrt(float64(nor2)))
		xyz.X *= invNor
		xyz.Y *= invNor
		xyz.Z *= invNor
	}
}

func (xyz XYZGenericType) Normalized() *XYZGenericType {
	nor2 := xyz.LengthSquared()
	if nor2 > 0 {
		invNor := GenericType(1) / GenericType(math.Sqrt(float64(nor2)))
		xyz.X *= invNor
		xyz.Y *= invNor
		xyz.Z *= invNor
	}
	return &xyz
}

func (xyz *XYZGenericType) Equals(other *XYZGenericType) bool {
	return xyz.X == other.X && xyz.Y == other.Y && xyz.Z == other.Z
}

func (xyz *XYZGenericType) NotEquals(other *XYZGenericType) bool {
	return xyz.X != other.X || xyz.Y != other.Y || xyz.Z != other.Z
}
