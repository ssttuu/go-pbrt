//go:generate mockgen -source=scene.go -destination=scene.mock.go -package=pbrt

package pbrt

type Scene interface {
	Aggregate() Aggregate
	Intersect(r *Ray, si *SurfaceInteraction) bool
	IntersectP(r *Ray) bool
	IntersectTr(r *Ray, si *SurfaceInteraction, sampler Sampler, transmittance Spectrum) bool
	Light(index int) Light
	Lights() []Light
	InfiniteLights() []Light
}

func NewScene(aggregate Aggregate, lights, infiniteLights []Light) Scene {
	return &scene{
		lights:         lights,
		infiniteLights: infiniteLights,
		aggregate:      aggregate,
		worldBound:     aggregate.WorldBound(),
	}
}

type scene struct {
	lights         []Light
	infiniteLights []Light

	aggregate  Aggregate
	worldBound *Bounds3
}

func (s *scene) Aggregate() Aggregate {
	return s.aggregate
}

func (s *scene) Intersect(r *Ray, si *SurfaceInteraction) bool {
	return s.aggregate.Intersect(r, si)
}

func (s *scene) IntersectP(r *Ray) bool {
	return s.aggregate.IntersectP(r)
}

func (s *scene) IntersectTr(r *Ray, si *SurfaceInteraction, sampler Sampler, transmittance Spectrum) bool {
	transmittance.SetAll(1.0)
	for {
		hitSurface := s.Intersect(r, si)

		// Accumulate beam transmittance for ray segment
		if r.Medium != nil {
			transmittance.MulAssign(r.Medium.Tr(r, sampler))
		}

		// Initialize next ray segment or terminate transmittance computation
		if !hitSurface {
			return false
		}
		if si.Primitive != nil && si.Primitive.GetMaterial() != nil {
			return true
		}
		r = si.SpawnRay(r.Direction)
	}
}

func (s *scene) Light(index int) Light {
	return s.lights[index]
}

func (s *scene) Lights() []Light {
	return s.lights
}

func (s *scene) InfiniteLights() []Light {
	return s.infiniteLights
}
