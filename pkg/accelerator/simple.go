package accelerator

import (
	"log"

	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

func NewSimpleAggregate(primitives []pbrt.Primitive) *Simple {
	sa := &Simple{
		primitives: primitives,
	}
	if len(primitives) == 0 {
		return sa
	}
	sa.bounds = *primitives[0].WorldBound()
	for _, p := range primitives[1:] {
		sa.bounds.Union(p.WorldBound())
	}
	return sa
}

type Simple struct {
	primitives []pbrt.Primitive
	bounds     pbrt.Bounds3
}

func (a *Simple) GetAreaLight() pbrt.AreaLighter {
	log.Panic("aggregate.GetAreaLight called; should have gone to Geometric Primitive")
	return nil
}

func (a *Simple) GetMaterial() pbrt.Material {
	log.Panic("aggregate.GetMaterial called; should have gone to GeometricPrimitive")
	return nil
}

func (a *Simple) ComputeScatteringFunctions(si *pbrt.SurfaceInteraction, mode pbrt.TransportMode, allowMultipleLobes bool) {
	log.Panic("aggregate.ComputeScatteringFunctions called; should have gone to GeometricPrimitive")
}

func (a *Simple) WorldBound() *pbrt.Bounds3 {
	return &a.bounds
}

func (a *Simple) Intersect(r *pbrt.Ray, si *pbrt.SurfaceInteraction) bool {
	intersects := false
	minDist := math.Infinity
	closestInteraction := pbrt.NewSurfaceInteraction()

	for i := range a.primitives {
		interaction := pbrt.NewSurfaceInteraction()
		doesIntersect := a.primitives[i].Intersect(r, interaction)
		if doesIntersect {
			intersects = true
			dist := interaction.Point.Distance(r.Origin)
			if dist < minDist {
				minDist = dist
				*closestInteraction = *interaction
			}
		}
	}

	if intersects {
		*si = *closestInteraction
	}

	return intersects
}
func (a *Simple) IntersectP(r *pbrt.Ray) bool {
	for _, p := range a.primitives {
		intersects := p.IntersectP(r)
		if intersects {
			return true
		}
	}
	return false
}
