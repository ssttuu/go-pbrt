package pbrt
//
//type XGetter interface {
//	GetX() float64
//}
//
//type YGetter interface {
//	GetY() float64
//}
//
//type ZGetter interface {
//	GetZ() float64
//}
//
//type XSetter interface {
//	SetX(value float64)
//}
//
//type YSetter interface {
//	SetY(value float64)
//}
//
//type ZSetter interface {
//	SetZ(value float64)
//}
//
//// XY
//
//type XYGetter interface {
//	XGetter
//	YGetter
//}
//
//type XYSetter interface {
//	XSetter
//	YSetter
//}
//
//type XYAdder interface {
//	Add(other XYGetter) XYer
//	AddAssign(other XYGetter)
//}
//
//type XYSubtracter interface {
//	Sub(other XYGetter) XYer
//	SubAssign(other XYGetter)
//}
//
//type XYMultiplier interface {
//	Mul(other XYGetter) XYer
//	MulAssign(other XYGetter)
//}
//
//type XYDivider interface {
//	Div(other XYGetter) XYer
//	DivAssign(other XYGetter)
//}
//
//type XYDotter interface {
//	Dot(other XYGetter) float64
//}
//
//type XYNormalizer interface {
//	Normalize()
//}
//
//type XYer interface {
//	XYGetter
//	XYSetter
//	XYAdder
//	XYSubtracter
//	XYMultiplier
//	XYDivider
//	XYDotter
//	XYNormalizer
//}
//
//// XYZ
//
//type XYZGetter interface {
//	XGetter
//	YGetter
//	ZGetter
//}
//
//type XYZSetter interface {
//	XSetter
//	YSetter
//	ZSetter
//}
//
//type XYZAdder interface {
//	Add(other XYZGetter) XYZer
//	AddAssign(other XYZGetter)
//}
//
//type XYZSubtracter interface {
//	Sub(other XYZGetter) XYZer
//	SubAssign(other XYZGetter)
//}
//
//type XYZMultiplier interface {
//	Mul(other XYZGetter) XYZer
//	MulAssign(other XYZGetter)
//}
//
//type XYZDivider interface {
//	Div(other XYZGetter) XYZer
//	DivAssign(other XYZGetter)
//}
//
//type XYZDotter interface {
//	Dot(other XYZGetter) float64
//}
//
//type XYZNormalizer interface {
//	Normalize()
//}
//
//type XYZer interface {
//	XYZGetter
//	XYZSetter
//	XYZAdder
//	XYZSubtracter
//	XYZMultiplier
//	XYZDivider
//}
