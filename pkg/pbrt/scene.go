package pbrt


type Scene struct {
	lights []Light
	infiniteLights []Light

	aggregate Primitive
	worldBound *Bounds3
}

func (s *Scene) Intersect(r *Ray, si *SurfaceInteraction) bool {
	return s.aggregate.Intersect(r, si)
}

func (s *Scene) IntersectP(r *Ray) bool {
	return s.aggregate.IntersectP(r)
}

func (s *Scene) IntersectTr(r *Ray, sampler *Sampler, si *SurfaceInteraction, transmittance *Spectrum) bool {

	*transmittance = NewSpectrum(1.0)
	for {
		hitSurface := s.Intersect(r, si)
		if r.medium != nil {
			transmittance.Mul(r.medium.Tr(r, sampler))
		}

		if !hitSurface {
			return false
		}
		if si.primitive.GetMaterial() != nil {
			return true
		}
		r = si.SpawnRay(r.direction)
	}
}