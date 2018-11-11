package pbrt

import (
	"log"
)

type Primitiver interface {
	WorldBound() *Bounds3
	Intersect(r *Ray) (bool, *SurfaceInteraction)
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() Materialer
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
}

type Aggregator interface {
	WorldBound() *Bounds3
	Intersect(r *Ray) (bool, *SurfaceInteraction)
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() Materialer
}

type GeometricPrimitive struct {
	Shape          Shaper
	material       Materialer
	areaLight      AreaLighter
	mediumAccessor *MediumAccessor
}

func NewGeometricPrimitive(shape Shaper, m Materialer) *GeometricPrimitive {
	return &GeometricPrimitive{
		Shape:          shape,
		material:       m,
		areaLight:      nil,
		mediumAccessor: new(MediumAccessor),
	}
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

func (p *GeometricPrimitive) GetAreaLight() AreaLighter {
	return p.areaLight
}

func (p *GeometricPrimitive) GetMaterial() Materialer {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (p *GeometricPrimitive) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	if p.material == nil {
		log.Panicf("NO MATERIAL: %+v", p.Shape.GetName())
	}

	p.material.ComputeScatteringFunctions(si, mode, allowMultipleLobes)
}

//func NewGeometricPrimitive(Shape *Shape, material Material, areaLigh)

type TransformedPrimitive struct {
	primitive        Primitiver
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
	log.Panic("Aggregate.GetAreaLight called; should have gone to Geometric Primitiver")
	return nil
}

func (a *Aggregate) GetMaterial() Materialer {
	log.Panic("Aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

//func (a *Aggregate) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
//	log.Panic("Aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
//}

type SimpleAggregate struct {
	*Aggregate

	primitives []Primitiver
	bounds     *Bounds3
}

func NewSimpleAggregate(primitives []Primitiver) *SimpleAggregate {
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
	intersects := false
	minDist := Infinity
	var minIsect *SurfaceInteraction

	for _, p := range a.primitives {
		doesIntersect, isect := p.Intersect(r)
		if doesIntersect {
			intersects = true
			dist2 := isect.point.Distance(r.origin)
			if dist2 < minDist {
				minDist = dist2
				minIsect = isect
			}
		}
	}

	return intersects, minIsect
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
