package pbrt

import (
	"fmt"
	"math"
)

type Normal3 struct {
	x float64
	y float64
	z float64
}

func NewNormal3(x, y, z float64) *Normal3 {
	return &Normal3{
		x: x,
		y: y,
		z: z,
	}
}

func (n *Normal3) ToVector3() *Vector3 {
	return &Vector3{
		x: n.x,
		y: n.y,
		z: n.z,
	}
}

func (n *Normal3) Clone() *Normal3 {
	return &Normal3{
		x: n.x,
		y: n.y,
		z: n.z,
	}
}

func (n *Normal3) GetX() float64 {
	return n.x
}

func (n *Normal3) GetY() float64 {
	return n.y
}

func (n *Normal3) GetZ() float64 {
	return n.z
}

func (n *Normal3) SetX(value float64) {
	n.x = value
}

func (n *Normal3) SetY(value float64) {
	n.y = value
}

func (n *Normal3) SetZ(value float64) {
	n.z = value
}

func (n *Normal3) Index(i int) float64 {
	switch i {
	case 0:
		return n.x
	case 1:
		return n.y
	default:
		return n.z
	}
}

func (n *Normal3) SetIndex(i int, v float64) {
	switch i {
	case 0:
		n.x = v
	case 1:
		n.y = v
	default:
		n.z = v
	}
}

func (n *Normal3) LengthSquared() float64 {
	return n.x * n.x + n.y * n.y + n.z * n.z
}

func (n *Normal3) Length() float64 {
	return math.Sqrt(n.LengthSquared())
}

func (n *Normal3) String() string {
	return fmt.Sprintf("Normal3(%f.3, %f.3, %f.3)", n.x, n.y, n.z)
}

func (n *Normal3) Add(other XYZGetter) Normal3 {
	return Normal3{
		x: n.x + other.GetX(),
		y: n.y + other.GetY(),
		z: n.z + other.GetZ(),
	}
}

func (n *Normal3) AddAssign(other XYZGetter) {
	n.x += other.GetX()
	n.y += other.GetY()
	n.z += other.GetZ()
}

func (n *Normal3) Sub(other XYZGetter) Normal3 {
	return Normal3{
		x: n.x - other.GetX(),
		y: n.y - other.GetY(),
		z: n.z - other.GetZ(),
	}
}

func (n *Normal3) SubAssign(other XYZGetter) {
	n.x += other.GetX()
	n.y += other.GetY()
	n.z += other.GetZ()
}

func (n *Normal3) Mul(other XYZGetter) Normal3 {
	return Normal3{
		x: n.x * other.GetX(),
		y: n.y * other.GetY(),
		z: n.z * other.GetZ(),
	}
}

func (n *Normal3) MulAssign(other XYZGetter) {
	n.x += other.GetX()
	n.y += other.GetY()
	n.z += other.GetZ()
}

func (n *Normal3) MulScalar(other float64) *Normal3 {
	return &Normal3{n.x * other, n.y * other, n.z * other}
}

func (n *Normal3) Div(other XYZGetter) Normal3 {
	return Normal3{
		x: n.x / other.GetX(),
		y: n.y / other.GetY(),
		z: n.z / other.GetZ(),
	}
}

func (n *Normal3) DivAssign(other XYZGetter) {
	n.x += other.GetX()
	n.y += other.GetY()
	n.z += other.GetZ()
}

func (n *Normal3) Dot(other XYZGetter) float64 {
	return n.GetX() * other.GetX() + n.GetY() * other.GetY() + n.GetZ() * other.GetZ()
}

func (n *Normal3) Normalize() {
	nor2 := n.LengthSquared()
	if (nor2 > 0) {
		invNor := 1 / math.Sqrt(nor2)
		n.x *= invNor
		n.y *= invNor
		n.z *= invNor
	}
}

func (n *Normal3) Abs() *Normal3 {
	return &Normal3{
		x: math.Abs(n.x),
		y: math.Abs(n.y),
		z: math.Abs(n.z),
	}
}

func (n *Normal3) Equals(other XYZGetter) bool {
	return n.x == other.GetX() && n.y == other.GetY() && n.z == other.GetZ()
}

func (n *Normal3) NotEquals(other XYZGetter) bool {
	return n.x != other.GetX() || n.y != other.GetY() || n.z != other.GetZ()
}
