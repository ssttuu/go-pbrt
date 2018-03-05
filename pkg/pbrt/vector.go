package pbrt

import (
	"math"
	"fmt"
)


type Vector2 struct {
	x float64
	y float64
}

func (v *Vector2) GetX() float64 {
	return v.x
}

func (v *Vector2) GetY() float64 {
	return v.y
}

func (v *Vector2) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y
}

func (v *Vector2) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *Vector2) String() string {
	return fmt.Sprintf("Vector2(%f.3, %f.3)", v.x, v.y)
}

func (v *Vector2) Add(other XYGetter) *Vector2 {
	return &Vector2{v.x + other.GetX(), v.y + other.GetY()}
}

func (v *Vector2) AddAssign(other XYGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
}

func (v *Vector2) Sub(other XYGetter) *Vector2 {
	return &Vector2{v.x - other.GetX(), v.y - other.GetY()}
}

func (v *Vector2) SubAssign(other XYGetter) {
	v.x -= other.GetX()
	v.y -= other.GetY()
}

func (v *Vector2) Mul(other XYGetter) *Vector2 {
	return &Vector2{v.x * other.GetX(), v.y * other.GetY()}
}

func (v *Vector2) MulAssign(other XYGetter) {
	v.x *= other.GetX()
	v.y *= other.GetY()
}

func (v *Vector2) MulScalar(other float64) *Vector2 {
	return &Vector2{v.x * other, v.y * other}
}

func (v *Vector2) Div(other XYGetter) *Vector2 {
	return &Vector2{v.x / other.GetX(), v.y / other.GetY()}
}

func (v *Vector2) DivAssign(other XYGetter) {
	v.x /= other.GetX()
	v.y /= other.GetY()
}

func (v *Vector2) Equals(other XYGetter) bool {
	return v.x == other.GetX() && v.y == other.GetY()
}

func (v *Vector2) NotEquals(other XYGetter) bool {
	return v.x != other.GetX() || v.y != other.GetY()
}

type Vector3 struct {
	x float64
	y float64
	z float64
}

func NewVector3(x, y, z float64) *Vector3 {
	return &Vector3{
		x: x,
		y: y,
		z: z,
	}
}

func (v *Vector3) Clone() *Vector3 {
	return &Vector3{
		x: v.x,
		y: v.y,
		z: v.z,
	}
}

func (v *Vector3) GetX() float64 {
	return v.x
}

func (v *Vector3) GetY() float64 {
	return v.y
}

func (v *Vector3) GetZ() float64 {
	return v.z
}

func (v *Vector3) SetX(value float64) {
	v.x = value
}

func (v *Vector3) SetY(value float64) {
	v.y = value
}

func (v *Vector3) SetZ(value float64) {
	v.z = value
}

func (v *Vector3) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y + v.z * v.z
}

func (v *Vector3) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *Vector3) String() string {
	return fmt.Sprintf("Vector3(%f.3, %f.3, %f.3)", v.x, v.y, v.z)
}

func (v *Vector3) Add(other XYZGetter) *Vector3 {
	return &Vector3{
		x: v.x + other.GetX(),
		y: v.y + other.GetY(),
		z: v.z + other.GetZ(),
	}
}

func (v *Vector3) AddAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Vector3) Sub(other XYZGetter) *Vector3 {
	return &Vector3{
		x: v.x - other.GetX(),
		y: v.y - other.GetY(),
		z: v.z - other.GetZ(),
	}
}

func (v *Vector3) SubAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Vector3) Mul(other XYZGetter) *Vector3 {
	return &Vector3{
		x: v.x * other.GetX(),
		y: v.y * other.GetY(),
		z: v.z * other.GetZ(),
	}
}

func (v *Vector3) MulAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Vector3) MulScalar(other float64) *Vector3 {
	return &Vector3{
		x: v.x * other,
		y: v.y * other,
		z: v.z * other,
	}
}

func (v *Vector3) Div(other XYZGetter) *Vector3 {
	return &Vector3{
		x: v.x / other.GetX(),
		y: v.y / other.GetY(),
		z: v.z / other.GetZ(),
	}
}

func (v *Vector3) DivAssign(other XYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *Vector3) Dot(other *Vector3) float64 {
	return v.x * other.x + v.y * other.y + v.z * other.z
}

func (v *Vector3) Cross(other *Vector3) *Vector3 {
	return &Vector3{(v.y * other.z) - (v.z * other.y), (v.z * other.x) - (v.x * other.z), (v.x * other.y) - (v.y * other.x)}
}

func (v *Vector3) Normalize() {
	nor2 := v.LengthSquared()
	if (nor2 > 0) {
		invNor := 1 / math.Sqrt(nor2)
		v.x *= invNor
		v.y *= invNor
		v.z *= invNor
	}
}

func (v Vector3) Normalized() *Vector3 {
	nor2 := v.LengthSquared()
	if (nor2 > 0) {
		invNor := 1 / math.Sqrt(nor2)
		v.x *= invNor
		v.y *= invNor
		v.z *= invNor
	}
	return &v
}

func (v *Vector3) Equals(other XYZGetter) bool {
	return v.x == other.GetX() && v.y == other.GetY() && v.z == other.GetZ()
}

func (v *Vector3) NotEquals(other XYZGetter) bool {
	return v.x != other.GetX() || v.y != other.GetY() || v.z != other.GetZ()
}

