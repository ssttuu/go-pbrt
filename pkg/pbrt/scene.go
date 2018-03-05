package pbrt


type Scene struct {
	lights []Light
	infiniteLights []Light

	aggregate Primitive
	worldBound *Bounds3
}

func (s *Scene) Intersect(r *Ray, si SurfaceInteraction) bool {
	return s.aggregate.Intersect(r, si)
}

func (s *Scene) IntersectP(r *Ray) bool {
	return s.aggregate.IntersectP(r)
}

//func (s *Scene) IntersectTr(r *Ray, sampler *Sampler, si *SurfaceInteraction, transmittance *Spectrum) bool {
//
//	*transmittance = Spectrum{}
//}