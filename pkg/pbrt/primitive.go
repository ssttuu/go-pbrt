//go:generate mockgen -source=primitive.go -destination=primitive.mock.go -package=pbrt

package pbrt

import (
	"log"
)

type Primitive interface {
	Intersect(r *Ray, si *SurfaceInteraction) bool
	IntersectP(r *Ray) bool
	GetAreaLight() AreaLighter
	GetMaterial() Material
	ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool)
	WorldBound() *Bounds3
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
	r.TMax = tHit
	si.Primitive = p

	if p.mediumAccessor.IsMediumTransition() {
		si.MediumAccessor = p.mediumAccessor
	} else {
		si.MediumAccessor = &MediumAccessor{r.Medium, r.Medium}
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
		log.Panicf("NO MATERIAL: %+v", p.Shape)
	}

	p.material.ComputeScatteringFunctions(si, mode, allowMultipleLobes)
}

//func NewGeometricPrimitive(shape *shape, material material, areaLigh)

func NewTransformedPrimitive(p Primitive, p2w *AnimatedTransform) *TransformedPrimitive {
	return &TransformedPrimitive{
		primitive:        p,
		primitiveToWorld: p2w,
	}
}

type TransformedPrimitive struct {
	primitive        Primitive
	primitiveToWorld *AnimatedTransform
}

func (p *TransformedPrimitive) Intersect(r *Ray, si *SurfaceInteraction) bool {
	interpolatedPrimToWorld := p.primitiveToWorld.Interpolate(r.Time)
	ray, _, _ := interpolatedPrimToWorld.Inverse().TransformRay(r)

	intersects := p.primitive.Intersect(ray, si)
	if !intersects {
		return false
	}
	r.TMax = ray.TMax

	if !interpolatedPrimToWorld.IsIdentity() {
		si = interpolatedPrimToWorld.TransformSurfaceInteraction(si)
	}

	return true
}

func (p *TransformedPrimitive) IntersectP(r *Ray) bool {
	interpolatedPrimToWorld := p.primitiveToWorld.Interpolate(r.Time)
	ray, _, _ := interpolatedPrimToWorld.Inverse().TransformRay(r)
	return p.primitive.IntersectP(ray)
}

func (p *TransformedPrimitive) GetAreaLight() AreaLighter {
	return nil
}
func (p *TransformedPrimitive) GetMaterial() Material {
	return nil
}
func (p *TransformedPrimitive) ComputeScatteringFunctions(si *SurfaceInteraction, mode TransportMode, allowMultipleLobes bool) {
	log.Panic("should not be called")
}

func (p *TransformedPrimitive) WorldBound() *Bounds3 {
	return p.primitiveToWorld.MotionBounds(p.primitive.WorldBound())
}
