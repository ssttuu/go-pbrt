package pbrt


type Scene struct {
	lights []Light
	infiniteLights []Light

	aggregate Primitive
	worldBound *Bounds3
}

func (s *Scene) Intersect(ray *Ray, si SurfaceInteraction) bool {
	return s.aggregate.Intersect(ray, si)
}

func (s *Scene) IntersectP(ray *Ray) bool {
	return s.aggregate.IntersectP(ray)
}
