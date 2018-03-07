package pbrt

import "log"

type Primitive interface {
	WorldBound() Bounds3
	Intersect(r *Ray, si *SurfaceInteraction) bool
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() *Material
	ComputeScatteringFunctions(si *SurfaceInteraction, arena *MemoryArena, mode TransportMode, allowMultipleLobes bool)
}

type GeometricPrimitive struct {
	shape          *Shape
	material       *Material
	areaLight      *AreaLight
	mediumAccessor *MediumAccessor
}

type TransformedPrimitive struct {
	primitive        Primitive
	primitiveToWorld AnimatedTransform
}

func (tp *TransformedPrimitive) Intersect(r *Ray, si *SurfaceInteraction) bool {
	interpolatedPrimToWorld := tp.primitiveToWorld.Interpolate(r.time)
	ray := interpolatedPrimToWorld.Inverse().TransformRay(r)
	if !tp.primitive.Intersect(ray, si) {
		return false
	}
	r.tMax = ray.tMax

	if !interpolatedPrimToWorld.IsIdentity() {
		si = interpolatedPrimToWorld.TransformSurfaceInteraction(si)
	}

	return true
}

func (tp *TransformedPrimitive) IntersectP(r *Ray) bool {
	interpolatedPrimToWorld := tp.primitiveToWorld.Interpolate(r.time).Inverse()
	return tp.primitive.IntersectP(interpolatedPrimToWorld.TransformRay(r))
}

type Aggregate struct {

}

func (a *Aggregate) GetAreaLight() *AreaLight {
	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitive")
	return nil
}

func (a *Aggregate) GetMaterial() *Material {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (a *Aggregate) ComputeScatteringFunctions(si *SurfaceInteraction, arena *MemoryArena, mode TransportMode, allowMultipleLobes bool) {
	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
}