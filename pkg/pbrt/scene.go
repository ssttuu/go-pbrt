package pbrt

type Scene struct {
	Lights         []Light
	infiniteLights []Light

	Aggregate  Aggregate
	WorldBound *Bounds3
}

func (s *Scene) Intersect(r *Ray, si *SurfaceInteraction) bool {
	return s.Aggregate.Intersect(r, si)
}

func (s *Scene) IntersectP(r *Ray) bool {
	return s.Aggregate.IntersectP(r)
}

func (s *Scene) IntersectTr(r *Ray, si *SurfaceInteraction, sampler Sampler, transmittance Spectrum) bool {
	transmittance.SetAll(1.0)
	for {
		hitSurface := s.Intersect(r, si)

		// Accumulate beam transmittance for ray segment
		if r.medium != nil {
			transmittance.MulAssign(r.medium.Tr(r, sampler))
		}

		// Initialize next ray segment or terminate transmittance computation
		if !hitSurface {
			return false
		}
		if si.Primitive != nil && si.Primitive.GetMaterial() != nil {
			return true
		}
		r = si.SpawnRay(r.direction)
	}
}
