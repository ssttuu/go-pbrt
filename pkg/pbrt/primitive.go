//go:generate mockgen -source=primitive.go -destination=primitive.mock.go -package=pbrt

package pbrt

import (
	"log"
)

type Primitive interface {
	WorldBound() *Bounds3
	Intersect(r *Ray, si *SurfaceInteraction) bool
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() Material
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
}

type Aggregate interface {
	Primitive
}

type GeometricPrimitive struct {
	Shape          Shape
	material       Material
	areaLight      AreaLighter
	mediumAccessor *MediumAccessor
}

func NewGeometricPrimitive(shape Shape, m Material) *GeometricPrimitive {
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

func (p *GeometricPrimitive) Intersect(r *Ray, si *SurfaceInteraction) bool {
	intersects, tHit := p.Shape.Intersect(r, si, true)
	if !intersects {
		return false
	}
	r.tMax = tHit
	si.Primitive = p

	if p.mediumAccessor.IsMediumTransition() {
		si.mediumAccessor = p.mediumAccessor
	} else {
		si.mediumAccessor = &MediumAccessor{r.medium, r.medium}
	}

	return true
}

func (p *GeometricPrimitive) GetAreaLight() AreaLighter {
	return p.areaLight
}

func (p *GeometricPrimitive) GetMaterial() Material {
	log.Panic("aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (p *GeometricPrimitive) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	if p.material == nil {
		log.Panicf("NO MATERIAL: %+v", p.Shape.GetName())
	}

	p.material.ComputeScatteringFunctions(si, mode, allowMultipleLobes)
}

//func NewGeometricPrimitive(shape *shape, material material, areaLigh)

type TransformedPrimitive struct {
	primitive        Primitive
	primitiveToWorld AnimatedTransform
}

func (tp *TransformedPrimitive) Intersect(r *Ray, si *SurfaceInteraction) bool {
	interpolatedPrimToWorld := tp.primitiveToWorld.Interpolate(r.time)
	ray := interpolatedPrimToWorld.Inverse().TransformRay(r)

	intersects := tp.primitive.Intersect(ray, si)
	if !intersects {
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
