package pbrt

import (
	"math"
	"fmt"
	"log"
)

type Integrator interface {
	GetSampler() Sampler
	GetCamera() Cameraer

	Preprocess(scene *Scene, sampler Sampler)
	//Render(scene *Scene)
	Li(ray *RayDifferential, scene *Scene, sampler Sampler, depth int) Spectrum
	SpecularReflect(ray *RayDifferential, si *SurfaceInteraction, scene *Scene, sampler Sampler, depth int) Spectrum
	SpecularTransmit(ray *RayDifferential, si *SurfaceInteraction, scene *Scene, sampler Sampler, depth int) Spectrum
}

func UniformSampleAllLights(it Interaction, scene *Scene, sampler Sampler, nLightSamples []int, handleMedia bool) Spectrum {
	L := NewSpectrum(0)
	for j, light := range scene.Lights {
		// accumulate contribution of jth light to L
		//light := scene.Lights[j]
		nSamples := nLightSamples[j]
		uLightSlice := sampler.Get2DArray(nSamples)
		uScatteringSlice := sampler.Get2DArray(nSamples)
		if len(uLightSlice) == 0 || len(uScatteringSlice) == 0 {
			fmt.Println("Zero length slice")
			// use a single sample for illumination from light
			uLight := sampler.Get2D()
			uScattering := sampler.Get2D()
			L.AddAssign(EstimateDirect(it, uScattering, light, uLight, scene, sampler, handleMedia, false))
		} else {
			Ld := NewSpectrum(0)
			for k := 0; k < nSamples; k++ {
				Ld.AddAssign(EstimateDirect(it, uScatteringSlice[k], light, uLightSlice[k], scene, sampler, handleMedia, false))
			}

			L.AddAssign(Ld.DivScalar(float64(nSamples)))
		}
	}
	return L
}

func UniformSampleOneLight(it Interaction, scene *Scene, sampler Sampler, handleMedia bool, lightDistrib *Distribution1D) Spectrum {
	// randomly choose a single light to sample
	nLights := len(scene.Lights)
	if nLights == 0 {
		return NewSpectrum(0)
	}
	var lightNum int
	var lightPdf float64
	if lightDistrib != nil {
		lightNum, lightPdf, _ = lightDistrib.SampleDiscrete(sampler.Get1D())
		if lightPdf == 0.0 {
			return NewSpectrum(0)
		}
	} else {
		lightNum = int(math.Min(sampler.Get1D()*float64(nLights), float64(nLights-1)))
		lightPdf = 1.0 / float64(nLights)
	}

	light := scene.Lights[lightNum]
	uLight := sampler.Get2D()
	uScattering := sampler.Get2D()

	spectrum := EstimateDirect(it, uScattering, light, uLight, scene, sampler, handleMedia, false)
	spectrum.DivScalar(lightPdf)
	return spectrum
}

func EstimateDirect(it Interaction, uScattering *Point2f, light Lighter, uLight *Point2f, scene *Scene, sampler Sampler, handleMedia bool, specular bool) Spectrum {
	bsdfFlags := BxDFType(BSDFAll & ^BSDFSpecular)
	if specular {
		bsdfFlags = BSDFAll
	}

	scatteringPdf := 0.0
	Ld := NewSpectrum(0)
	fmt.Printf("%+v\n", uLight)
	Li, wi, lightPdf, visibility := light.SampleLi(it, uLight)
	if lightPdf > 0 || !Li.IsBlack() {
		fmt.Printf("%+v, %+v, %+v, %+v \n", Li, wi, lightPdf, visibility)
	}
	return Li
	if lightPdf > 0 && !Li.IsBlack() {
		// compute BSDF for light sampling strategy
		var f Spectrum

		switch intr := it.(type) {
		case *SurfaceInteraction:
			fmt.Printf("SurfaceInteraction: %+v %+v\n", intr.wo, wi)
			f = intr.bsdf.F(intr.wo, wi, bsdfFlags).MulScalar(wi.Dot(intr.shading.normal))
			scatteringPdf = intr.bsdf.Pdf(intr.wo, wi, bsdfFlags)
		case *MediumInteraction:
			fmt.Printf("MediumInteraction: %+v\n", intr)
			p := intr.phase.P(intr.wo, wi)
			f = NewSpectrum(p)
			scatteringPdf = p
		default:
			// TODO: return error
			log.Panicf("UnknownInteraction: %+v\n", intr)
		}

		if !f.IsBlack() {
			// compute effect of visibility for light source sample
			if handleMedia {
				Li.MulAssign(visibility.Tr(scene, sampler))
			} else {
				if !visibility.Unoccluded(scene) {
					// Shadow ray blocked
					Li = NewSpectrum(0)
				} else {
					// Shadow ray unoccluded
				}
			}

			// add light's contribution to reflected radiance
			if !Li.IsBlack() {
				if IsDeltaLight(light.GetFlags()) {
					Ld.AddAssign(f.Mul(Li).DivScalar(lightPdf))
				} else {
					weight := PowerHeuristic(1.0, lightPdf, 1.0, scatteringPdf)
					Ld.AddAssign(f.Mul(Li).MulScalar(weight).DivScalar(lightPdf))
				}
			}
		}

		// sample SDF with multiple importance sampling
		if !IsDeltaLight(light.GetFlags()) {
			var f Spectrum
			sampledSpecular := false
			switch intr := it.(type) {
			case SurfaceInteraction:
				// sample scattered direction for surface interactions
				var sampledType BxDFType
				f, wi, scatteringPdf, sampledType = intr.bsdf.SampleF(intr.wo, uScattering, bsdfFlags)
				f.MulScalar(wi.AbsDot(intr.shading.normal))
				sampledSpecular = sampledType&BSDFSpecular != 0
			case MediumInteraction:
				// sample scattered direction for medium interactions
				p := intr.phase.SampleP(intr.wo, wi, uScattering)
				f = NewSpectrum(p)
				scatteringPdf = p
			}

			if !f.IsBlack() && scatteringPdf > 0.0 {
				// Account for light contributions along sampled direction wi
				weight := 1.0
				if !sampledSpecular {
					lightPdf, wi = light.PdfLi(it)
					if lightPdf == 0 {
						return Ld
					}
					weight = PowerHeuristic(1.0, scatteringPdf, 1.0, lightPdf)
				}

				// find intersection and compute transmittance
				ray := it.SpawnRay(wi)
				Tr := NewSpectrum(1.0)
				var foundSurfaceInteraction bool
				var lightIsect *SurfaceInteraction
				if handleMedia {
					foundSurfaceInteraction, lightIsect, Tr = scene.IntersectTr(ray, sampler)
				} else {
					foundSurfaceInteraction, lightIsect = scene.Intersect(ray)
				}

				// add light contribution from material sampling
				Li := NewSpectrum(0)
				if foundSurfaceInteraction {
					if lightIsect.primitive.GetAreaLight() == light {
						Li = lightIsect.Le(wi.MulScalar(-1))
					}
				} else {
					Li = light.Le(NewRayDifferentialFromRay(ray))
				}
				if !Li.IsBlack() {
					Ld.AddAssign(f.Mul(Li).Mul(Tr).MulScalar(weight / scatteringPdf))
				}
			}
		}
	}

	return Ld
}

type SamplerIntegrator struct {
	camera      Cameraer
	sampler     Sampler
	pixelBounds *Bounds2i
}

func NewSamplerIntegrator(camera Cameraer, sampler Sampler, pixelBounds *Bounds2i) *SamplerIntegrator {
	return &SamplerIntegrator{
		camera:      camera,
		sampler:     sampler,
		pixelBounds: pixelBounds,
	}
}

func (s *SamplerIntegrator) GetSampler() Sampler {
	return s.sampler
}

func (s *SamplerIntegrator) GetCamera() Cameraer {
	return s.camera
}

func (s *SamplerIntegrator) Preprocess(scene *Scene, sampler Sampler) {

}

func Render(s Integrator, scene *Scene) {
	s.Preprocess(scene, s.GetSampler())
	// render image tiles in parallel
	camera := s.GetCamera()
	film := camera.GetFilm()

	// compute number of tiles to use for parallel computing
	sampleBounds := camera.GetFilm().GetSampleBounds()
	sampleExtent := sampleBounds.Diagonal()
	var tileSize int64 = 16
	nTiles := &Point2i{(sampleExtent.X + tileSize - 1) / tileSize, (sampleExtent.Y + tileSize - 1) / tileSize}

	var tileY, tileX int64
	for tileY = 0; tileY < nTiles.Y; tileY++ {
		for tileX = 0; tileX < nTiles.X; tileX++ {
			tile := &Point2i{tileX, tileY}

			// get sampler instance for tile
			seed := uint64(tile.Y*nTiles.X + tile.X)
			tileSampler := s.GetSampler().Clone(seed)

			// compute sample bounds for tile
			//fmt.Println(sampleBounds, tileSize)
			x0 := sampleBounds.Min.X + tile.X*tileSize
			x1 := int64(math.Min(float64(x0+tileSize), float64(sampleBounds.Max.X)))
			y0 := sampleBounds.Min.Y + tile.Y*tileSize
			y1 := int64(math.Min(float64(y0+tileSize), float64(sampleBounds.Max.Y)))
			tileBounds := &Bounds2i{&Point2i{x0, y0}, &Point2i{x1, y1}}

			// get film tile
			//filmTile := camera.GetFilm().GetFilmTile(tileBounds)

			for pixelY := tileBounds.Min.Y; pixelY < tileBounds.Max.Y; pixelY++ {
				for pixelX := tileBounds.Min.X; pixelX < tileBounds.Max.X; pixelX++ {
					pixelPos := &Point2i{pixelX, pixelY}
					tileSampler.StartPixel(pixelPos)

					//for tileSampler.StartNextSample() {
					// initialize CameraSample for current sample
					cameraSample := tileSampler.GetCameraSample(pixelPos)

					// generate camera ray for current sample

					rayWeight, rd := camera.GenerateRayDifferential(cameraSample)
					rd.ScaleDifferentials(1.0 / math.Sqrt(float64(tileSampler.GetSamplesPerPixel())))
					//fmt.Printf("%+v %+v %+v\n", rayWeight, rd.origin, rd.direction)

					// evaluate radiance along camera ray
					L := NewSpectrum(0)
					if rayWeight > 0 {
						L = s.Li(rd, scene, tileSampler, 0)
					}

					// TODO: add errors and return them
					if L.HasNaNs() {
						L = NewSpectrum(0.1)
					} else if L.Y() < -1e-5 {
						L = NewSpectrum(0.5)
					} else if math.IsInf(L.Y(), 1) {
						L = NewSpectrum(1.0)
					}

					// add camera ray's contribution to image
					//filmTile.AddSample(cameraSample.pFilm, L, rayWeight)
					px := film.getPixel(pixelPos)
					px.value = [3]float64{
						L.Index(0),
						L.Index(1),
						L.Index(2),
					}

					if !L.IsBlack() {
						fmt.Printf("Not black: %+v\n", L)
						//fmt.Printf("%+v %+v %+v\n", pixelPos, cameraSample.pFilm, film.getPixel(pixelPos))
					}
					//}

					// merge image tile into film
					//camera.GetFilm().MergeFilmTile(filmTile)
				}
			}
		}
	}

	camera.GetFilm().WriteImage(1.0)
}

func (s *SamplerIntegrator) Li(ray *RayDifferential, scene *Scene, sampler Sampler, depth int) Spectrum {
	return nil
}

func (s *SamplerIntegrator) SpecularReflect(ray *RayDifferential, si *SurfaceInteraction, scene *Scene, sampler Sampler, depth int) Spectrum {
	wo := si.wo
	t := BxDFType(BSDFReflection | BSDFSpecular)
	f, wi, pdf, t := si.bsdf.SampleF(wo, sampler.Get2D(), t)

	// return contribution of specular reflection
	ns := si.shading.normal
	if pdf > 0 && !f.IsBlack() && wi.AbsDot(ns) != 0.0 {
		// compute ray differential rd for specular reflection
		rd := NewRayDifferentialFromRay(si.SpawnRay(wi))
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = si.point.Add(si.dpdx)
			rd.ryOrigin = si.point.Add(si.dpdy)

			// compute differential reflected directions
			dndx := si.shading.dndu.MulScalar(si.dudx).Add(si.shading.dndv.MulScalar(si.dvdx))
			dndy := si.shading.dndu.MulScalar(si.dudy).Add(si.shading.dndv.MulScalar(si.dvdy))
			dwodx := ray.rxDirection.MulScalar(-1).Sub(wo)
			dwody := ray.ryDirection.MulScalar(-1).Sub(wo)
			dDNdx := dwodx.Dot(ns) + wo.Dot(dndx)
			dDNdy := dwody.Dot(ns) + wo.Dot(dndy)
			rd.rxDirection = wi.Sub(dwodx).Add(dndx.MulScalar(wo.Dot(ns)).Add(ns.MulScalar(dDNdx)).MulScalar(2.0))
			rd.ryDirection = wi.Sub(dwody).Add(dndy.MulScalar(wo.Dot(ns)).Add(ns.MulScalar(dDNdy)).MulScalar(2.0))
		}
		return f.Mul(s.Li(rd, scene, sampler, depth+1)).MulScalar(wi.AbsDot(ns) / pdf)
	}

	return NewSpectrum(0)
}

func (s *SamplerIntegrator) SpecularTransmit(ray *RayDifferential, si *SurfaceInteraction, scene *Scene, sampler Sampler, depth int) Spectrum {
	wo := si.wo
	f, wi, pdf, _ := si.bsdf.SampleF(wo, sampler.Get2D(), BxDFType(BSDFTransmission|BSDFSpecular))

	p := si.point
	ns := si.shading.normal

	L := NewSpectrum(0)
	if pdf > 0 && !f.IsBlack() && wi.AbsDot(ns) != 0.0 {
		rd := NewRayDifferentialFromRay(si.SpawnRay(wi))
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = p.Add(si.dpdx)
			rd.ryOrigin = p.Add(si.dpdy)

			eta := si.bsdf.Eta
			w := wo.MulScalar(-1)
			if wo.Dot(ns) < 0 {
				eta = 1.0 / eta
			}

			dndx := si.shading.dndu.MulScalar(si.dudx).Add(si.shading.dndv.MulScalar(si.dvdx))
			dndy := si.shading.dndu.MulScalar(si.dudy).Add(si.shading.dndv.MulScalar(si.dvdy))
			dwodx := ray.rxDirection.MulScalar(-1).Sub(wo)
			dwody := ray.ryDirection.MulScalar(-1).Sub(wo)
			dDNdx := dwodx.Dot(ns) + wo.Dot(dndx)
			dDNdy := dwody.Dot(ns) + wo.Dot(dndy)

			mu := eta*w.Dot(ns) - wi.Dot(ns)
			dmudx := (eta - (eta*eta*w.Dot(ns))/wi.Dot(ns)) * dDNdx
			dmudy := (eta - (eta*eta*w.Dot(ns))/wi.Dot(ns)) * dDNdy

			rd.rxDirection = wi.Add(dwodx.MulScalar(eta)).Sub(dndx.MulScalar(mu).Add(ns.MulScalar(dmudx)))
			rd.ryDirection = wi.Add(dwody.MulScalar(eta)).Sub(dndy.MulScalar(mu).Add(ns.MulScalar(dmudy)))
		}
		L = f.Mul(s.Li(rd, scene, sampler, depth+1)).MulScalar(wi.AbsDot(ns) / pdf)
	}

	return L
}

// TODO: move DirectLightingIntegrator out of core pbrt

type LightStrategy uint8

const (
	UniformSampleAll LightStrategy = iota + 1
	UniformSampleOne
)

type DirectLightingIntegrator struct {
	*SamplerIntegrator

	strategy      LightStrategy
	maxDepth      int
	nLightSamples []int
}

func NewDirectLightingIntegrator(strategy LightStrategy, maxDepth int, camera Cameraer, sampler Sampler, pixelBounds *Bounds2i) *DirectLightingIntegrator {
	return &DirectLightingIntegrator{
		SamplerIntegrator: NewSamplerIntegrator(camera, sampler, pixelBounds),
		strategy:          strategy,
		maxDepth:          maxDepth,
	}
}

func (dli *DirectLightingIntegrator) Preprocess(scene *Scene, sampler Sampler) {
	if dli.strategy == UniformSampleAll {
		// compute number of samples to use for each light
		for i := 0; i < len(scene.Lights); i++ {
			dli.nLightSamples = append(dli.nLightSamples, sampler.RoundCount(scene.Lights[i].GetSamples()))
		}

		// request samples for sampling all Lights
		for i := 0; i < dli.maxDepth; i++ {
			for j := 0; j < len(scene.Lights); j++ {
				sampler.Request2DArray(dli.nLightSamples[j])
				sampler.Request2DArray(dli.nLightSamples[j])
			}
		}
	}
}

func (dli *DirectLightingIntegrator) Li(ray *RayDifferential, scene *Scene, sampler Sampler, depth int) Spectrum {
	L := NewSpectrum(0)

	// find closest ray intersection or return background radiance
	//fmt.Println("Ray Origin, Direction:", ray.Ray.origin, ray.Ray.direction)
	intersects, si := scene.Intersect(ray.Ray)
	if !intersects {
		//fmt.Println("Does not intersect")
		for i := 0; i < len(scene.Lights); i++ {
			L.AddAssign(scene.Lights[i].Le(NewRayDifferentialFromRay(ray.Ray)))
		}
		return L
	}

	//fmt.Printf("Intersects with %+v | %+v\n", si.shape.WorldBound(), si.point)

	// compute scattering function for surface interaction
	// TODO
	si.ComputeScatteringFunctions(ray, false, Radiance)
	if si.bsdf == nil {
		return dli.Li(NewRayDifferentialFromRay(si.SpawnRay(ray.direction)), scene, sampler, depth)
	}

	//compute emitted light if ray hit an area light source
	L.AddAssign(si.Le(si.wo))

	if len(scene.Lights) > 0 {
		// compute direct lighting for DirectLightingIntegrator integrator
		if dli.strategy == UniformSampleAll {
			L.AddAssign(UniformSampleAllLights(si, scene, sampler, dli.nLightSamples, false))
		} else {
			L.AddAssign(UniformSampleOneLight(si, scene, sampler, false, nil))
		}
	}
	if depth+1 < dli.maxDepth {
		// trace rays for specular reflection and refraction
		L.AddAssign(dli.SpecularReflect(ray, si, scene, sampler, depth+1))
		L.AddAssign(dli.SpecularTransmit(ray, si, scene, sampler, depth+1))
	}

	return L
}
