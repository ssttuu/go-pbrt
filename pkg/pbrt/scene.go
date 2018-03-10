package pbrt


type Scene struct {
	Lights         []Lighter
	infiniteLights []Lighter

	Aggregate  Primitive
	WorldBound *Bounds3
}


func (s *Scene) Intersect(r *Ray) (bool, *SurfaceInteraction) {
	return s.Aggregate.Intersect(r)
}

func (s *Scene) IntersectP(r *Ray) bool {
	return s.Aggregate.IntersectP(r)
}

func (s *Scene) IntersectTr(r *Ray, sampler Sampler) (intersects bool, si *SurfaceInteraction, transmittance Spectrum) {
	transmittance = NewRGBSpectrum(1.0, 1.0, 1.0)
	for {
		hitSurface, si := s.Intersect(r)
		if r.medium != nil {
			transmittance.MulAssign(r.medium.Tr(r, sampler))
		}

		if !hitSurface {
			return false, si, transmittance
		}
		if si.primitive.GetMaterial() != nil {
			return true, si, transmittance
		}
		r = si.SpawnRay(r.direction)
	}
}