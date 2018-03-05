package pbrt

type Primitive interface {
	WorldBound() Bounds3
	Intersect(ray *Ray, si SurfaceInteraction) bool
	IntersectP(ray *Ray) bool
}