package pbrt

import (
	"log"
)

type Primitive interface {
	WorldBound() *Bounds3
	Intersect(r *Ray) (bool, *SurfaceInteraction)
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() *Material
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
}

type GeometricPrimitive struct {
	Shape          Shaper
	material       *Material
	areaLight      *AreaLight
	mediumAccessor *MediumAccessor
}

func (p *GeometricPrimitive) WorldBound() *Bounds3 {
	return p.Shape.WorldBound()
}

func (p *GeometricPrimitive) IntersectP(r *Ray) bool {
	return p.Shape.IntersectP(r, true)
}

func (p *GeometricPrimitive) Intersect(r *Ray) (bool, *SurfaceInteraction) {
	intersects, tHit, isect := p.Shape.Intersect(r, true)
	if !intersects {
		return false, nil
	}
	r.tMax = tHit
	isect.primitive = p

	if p.mediumAccessor.IsMediumTransition() {
		isect.mediumAccessor = p.mediumAccessor
	} else {
		//isect.mediumAccessor = &MediumAccessor{r.medium}
	}

	return true, isect
}

func (a *GeometricPrimitive) GetAreaLight() AreaLighter {
	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitive")
	return nil
}

func (a *GeometricPrimitive) GetMaterial() *Material {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (a *GeometricPrimitive) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
}

//func NewGeometricPrimitive(Shape *Shape, material Material, areaLigh)

type TransformedPrimitive struct {
	primitive        Primitive
	primitiveToWorld AnimatedTransform
}

func (tp *TransformedPrimitive) Intersect(r *Ray) (bool, *SurfaceInteraction) {
	interpolatedPrimToWorld := tp.primitiveToWorld.Interpolate(r.time)
	ray := interpolatedPrimToWorld.Inverse().TransformRay(r)
	intersects, si := tp.primitive.Intersect(ray)
	if !intersects {
		return false, si
	}
	r.tMax = ray.tMax

	if !interpolatedPrimToWorld.IsIdentity() {
		si = interpolatedPrimToWorld.TransformSurfaceInteraction(si)
	}

	return true, si
}

func (tp *TransformedPrimitive) IntersectP(r *Ray) bool {
	interpolatedPrimToWorld := tp.primitiveToWorld.Interpolate(r.time).Inverse()
	return tp.primitive.IntersectP(interpolatedPrimToWorld.TransformRay(r))
}

type Aggregate struct {
}

func (a *Aggregate) WorldBound() *Bounds3 {
	return nil
}
func (a *Aggregate) Intersect(r *Ray) (bool, *SurfaceInteraction) {
	return false, nil
}
func (a *Aggregate) IntersectP(r *Ray) bool {
	return false
}

func (a *Aggregate) GetAreaLight() AreaLighter {
	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitive")
	return nil
}

func (a *Aggregate) GetMaterial() *Material {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (a *Aggregate) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
}

type SimpleAggregate struct {
	*Aggregate

	primitives []Primitive
	bounds     *Bounds3
}

func NewSimpleAggregate(primitives []Primitive) *SimpleAggregate {
	sa := &SimpleAggregate{
		primitives: primitives,
	}
	if len(primitives) == 0 {
		return sa
	}
	sa.bounds = primitives[0].WorldBound()
	for _, p := range primitives[1:] {
		sa.bounds.Union(p.WorldBound())
	}
	return sa
}

func (a *SimpleAggregate) WorldBound() *Bounds3 {
	return a.bounds
}

func (a *SimpleAggregate) Intersect(r *Ray) (bool, *SurfaceInteraction) {
	for _, p := range a.primitives {
		intersects, isect := p.Intersect(r)
		if intersects {
			return intersects, isect
		}
	}
	return false, nil
}
func (a *SimpleAggregate) IntersectP(r *Ray) bool {
	for _, p := range a.primitives {
		intersects := p.IntersectP(r)
		if intersects {
			return intersects
		}
	}
	return false
}
