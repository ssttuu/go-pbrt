package geometry

import (
	"math"
	"fmt"
)

type GetXer interface {
	GetX() float64
}

type GetYer interface {
	GetY() float64
}

type GetZer interface {
	GetZ() float64
}

type SetXer interface {
	SetX(value float64)
}

type SetYer interface {
	SetY(value float64)
}

type SetZer interface {
	SetZ(value float64)
}

type VectorXYZGetter interface {
	GetXer
	GetYer
	GetZer
}

type VectorXYZSetter interface {
	SetXer
	SetYer
	SetZer
}

type VectorXYZAdder interface {
	Add(other VectorXYZGetter)
}

type VectorXYZSubtracter interface {
	Subtract(other *VectorXYZ)
}

type VectorXYZMultiplier interface {
	Multiply(other VectorXYZGetter)
}

type VectorXYZDivider interface {
	Divide(other VectorXYZGetter)
}

type VectorXYZDotter interface {
	Dot(other *VectorXYZ) float64
}

type VectorXYZNormalizer interface {
	Normalize()
}

type VectorXYZCloner interface {
	Clone() VectorXYZer
}

type VectorXYer interface {
	GetXer
	GetYer
}

type VectorXYZer interface {
	VectorXYZCloner
	VectorXYZGetter
	VectorXYZSetter
	VectorXYZAdder
	VectorXYZSubtracter
	VectorXYZMultiplier
	VectorXYZDivider
	VectorXYZDotter
	VectorXYZNormalizer
}

func VectorXYZAdd(v1 VectorXYZCloner, v2 VectorXYZGetter) VectorXYZer {
	v3 := v1.Clone()
	v3.Add(v2)
	return v3
}

func VectorXYZSubtract(v1, v2 *VectorXYZ) *VectorXYZ {
	v3 := v1.Clone()
	v3.Subtract(v2)
	return v3
}

func VectorXYZMultiply(v1 VectorXYZCloner, v2 VectorXYZGetter) VectorXYZer {
	v3 := v1.Clone()
	v3.Multiply(v2)
	return v3
}

func VectorXYZDivide(v1 VectorXYZCloner, v2 VectorXYZGetter) VectorXYZer {
	v3 := v1.Clone()
	v3.Divide(v2)
	return v3
}

func VectorXYZEquals(v1, v2 VectorXYZGetter) bool {
	return v1.GetX() == v2.GetX() && v1.GetY() == v2.GetY() && v1.GetZ() == v2.GetZ()
}

func VectorXYZNotEquals(v1, v2 VectorXYZGetter) bool {
	return v1.GetX() != v2.GetX() || v1.GetY() != v2.GetY() || v1.GetZ() != v2.GetZ()
}

type VectorXY struct {
	x float64
	y float64
}

func (v *VectorXY) GetX() float64 {
	return v.x
}

func (v *VectorXY) GetY() float64 {
	return v.y
}

func (v *VectorXY) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y
}

func (v *VectorXY) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *VectorXY) String() string {
	return fmt.Sprintf("Vector2(%f.3, %f.3)", v.x, v.y)
}

func (v *VectorXY) Add(other VectorXYer) VectorXY {
	return VectorXY{v.x + other.GetX(), v.y + other.GetY()}
}

func (v *VectorXY) Subtract(other VectorXYer) VectorXY {
	return VectorXY{v.x - other.GetX(), v.y - other.GetY()}
}

func (v *VectorXY) Multiply(other VectorXYer) VectorXY {
	return VectorXY{v.x * other.GetX(), v.y * other.GetY()}
}

func (v *VectorXY) Divide(other VectorXYer) VectorXY {
	return VectorXY{v.x / other.GetX(), v.y / other.GetY()}
}

func (v *VectorXY) Equals(other VectorXYer) bool {
	return v.x == other.GetX() && v.y == other.GetY()
}

func (v *VectorXY) NotEquals(other VectorXYer) bool {
	return v.x != other.GetX() || v.y != other.GetY()
}

type VectorXYZ struct {
	x float64
	y float64
	z float64
}

func NewVectorXYZ(x, y, z float64) *VectorXYZ {
	return &VectorXYZ{
		x: x,
		y: y,
		z: z,
	}
}

func (v *VectorXYZ) Clone() *VectorXYZ {
	return &VectorXYZ{
		x: v.x,
		y: v.y,
		z: v.z,
	}
}

func (v *VectorXYZ) GetX() float64 {
	return v.x
}

func (v *VectorXYZ) GetY() float64 {
	return v.y
}

func (v *VectorXYZ) GetZ() float64 {
	return v.z
}

func (v *VectorXYZ) SetX(value float64) {
	v.x = value
}

func (v *VectorXYZ) SetY(value float64) {
	v.y = value
}

func (v *VectorXYZ) SetZ(value float64) {
	v.z = value
}

func (v *VectorXYZ) LengthSquared() float64 {
	return v.x * v.x + v.y * v.y + v.z * v.z
}

func (v *VectorXYZ) Length() float64 {
	return math.Sqrt(v.LengthSquared())
}

func (v *VectorXYZ) String() string {
	return fmt.Sprintf("Vector3(%f.3, %f.3, %f.3)", v.x, v.y, v.z)
}

func (v *VectorXYZ) Add(other VectorXYZGetter) {
	v.x += other.GetX()
	v.y += other.GetY()
	v.z += other.GetZ()
}

func (v *VectorXYZ) Subtract(other *VectorXYZ) {
	v.x -= other.x
	v.y -= other.y
	v.z -= other.z
}

func (v *VectorXYZ) Multiply(other VectorXYZGetter) {
	v.x *= other.GetX()
	v.y *= other.GetY()
	v.z *= other.GetZ()
}

func (v *VectorXYZ) Divide(other VectorXYZGetter) {
	v.x /= other.GetX()
	v.y /= other.GetY()
	v.z /= other.GetZ()
}

func (v *VectorXYZ) Dot(other *VectorXYZ) float64 {
	return v.x * other.x + v.y * other.y + v.z * other.z
}

func (v *VectorXYZ) Normalize() {
	nor2 := v.LengthSquared()
	if (nor2 > 0) {
		invNor := 1 / math.Sqrt(nor2)
		v.x *= invNor
		v.y *= invNor
		v.z *= invNor
	}
}

func (v *VectorXYZ) Equals(other VectorXYZGetter) bool {
	return v.x == other.GetX() && v.y == other.GetY() && v.z == other.GetZ()
}

func (v *VectorXYZ) NotEquals(other VectorXYZGetter) bool {
	return v.x != other.GetX() || v.y != other.GetY() || v.z != other.GetZ()
}

