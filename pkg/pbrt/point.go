package pbrt

import (
	"math"
	"fmt"
)


type Point2 struct {
	x float64
	y float64
}

func (v *Point2) GetX() float64 {
	return v.x
}

func (v *Point2) GetY() float64 {
	return v.y
}

func (v *Point2) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y
}

func (v *Point2) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *Point2) String() string {
	return fmt.Sprintf("Point2(%f.3, %f.3)", v.x, v.y)
}

func (v *Point2) Add(other XYGetter) *Point2 {
	return &Point2{v.x + other.GetX(), v.y + other.GetY()}
}

func (v *Point2) AddAssign(other XYGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
}

func (v *Point2) Sub(other XYGetter) *Vector2 {
	return &Vector2{v.x - other.GetX(), v.y - other.GetY()}
}

func (v *Point2) SubAssign(other XYGetter) {
	v.x -= other.GetX()
	v.y -= other.GetY()
}

func (v *Point2) Mul(other XYGetter) *Point2 {
	return &Point2{v.x * other.GetX(), v.y * other.GetY()}
}

func (v *Point2) MulAssign(other XYGetter) {
	v.x *= other.GetX()
	v.y *= other.GetY()
}

func (v *Point2) Div(other XYGetter) *Point2 {
	return &Point2{v.x / other.GetX(), v.y / other.GetY()}
}

func (v *Point2) DivAssign(other XYGetter) {
	v.x /= other.GetX()
	v.y /= other.GetY()
}

func (v *Point2) Equals(other XYGetter) bool {
	return v.x == other.GetX() && v.y == other.GetY()
}

func (v *Point2) NotEquals(other XYGetter) bool {
	return v.x != other.GetX() || v.y != other.GetY()
}

type Point3 struct {
	x float64
	y float64
	z float64
}

func NewPoint3(x, y, z float64) *Point3 {
	return &Point3{
		x: x,
		y: y,
		z: z,
	}
}

func (v *Point3) Clone() *Point3 {
	return &Point3{
		x: v.x,
		y: v.y,
		z: v.z,
	}
}

func (v *Point3) GetX() float64 {
	return v.x
}

func (v *Point3) GetY() float64 {
	return v.y
}

func (v *Point3) GetZ() float64 {
	return v.z
}

func (v *Point3) SetX(value float64) {
	v.x = value
}

func (v *Point3) SetY(value float64) {
	v.y = value
}

func (v *Point3) SetZ(value float64) {
	v.z = value
}

func (p *Point3) Index(i int) float64 {
	switch i {
	case 0:
		return p.x
	case 1:
		return p.y
	default:
		return p.z
	}
}

func (p *Point3) SetIndex(i int, v float64) {
	switch i {
	case 0:
		p.x = v
	case 1:
		p.y = v
	default:
		p.z = v
	}
}

func (v *Point3) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y + v.z * v.z
}

func (v *Point3) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *Point3) String() string {
	return fmt.Sprintf("Point3(%f.3, %f.3, %f.3)", v.x, v.y, v.z)
}

func (v *Point3) Add(other XYZGetter) *Point3 {
	return &Point3{
		x: v.x + other.GetX(),
		y: v.y + other.GetY(),
		z: v.z + other.GetZ(),
	}
}

func (v *Point3) AddAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Point3) Sub(other XYZGetter) *Vector3 {
	return &Vector3{
		x: v.x - other.GetX(),
		y: v.y - other.GetY(),
		z: v.z - other.GetZ(),
	}
}

func (v *Point3) SubAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Point3) Mul(other XYZGetter) *Point3 {
	return &Point3{
		x: v.x * other.GetX(),
		y: v.y * other.GetY(),
		z: v.z * other.GetZ(),
	}
}

func (v *Point3) MulAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Point3) Div(other XYZGetter) *Point3 {
	return &Point3{
		x: v.x / other.GetX(),
		y: v.y / other.GetY(),
		z: v.z / other.GetZ(),
	}
}

func (v *Point3) DivAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Point3) Dot(other XYZGetter) float64 {
	return v.GetX() * other.GetX() + v.GetY() * other.GetY() + v.GetZ() * other.GetZ()
}

func (v *Point3) Normalize() {
	nor2 := v.LengthSquared()
	if (nor2 > 0) {
		invNor := 1 / math.Sqrt(nor2)
		v.x *= invNor
		v.y *= invNor
		v.z *= invNor
	}
}

func (v *Point3) Equals(other XYZGetter) bool {
	return v.x == other.GetX() && v.y == other.GetY() && v.z == other.GetZ()
}

func (v *Point3) NotEquals(other XYZGetter) bool {
	return v.x != other.GetX() || v.y != other.GetY() || v.z != other.GetZ()
}

