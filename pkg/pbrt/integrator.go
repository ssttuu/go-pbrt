package pbrt

import (
	"context"
	"log"
	"math"

	"github.com/pkg/errors"
	"golang.org/x/sync/errgroup"
)

type Integrator interface {
	GetSampler() Sampler
	GetCamera() Camera

	Preprocess(scene Scene, sampler Sampler)
	//Render(scene Scene)
	Li(ctx context.Context, ray *RayDifferential, scene Scene, sampler Sampler, depth int) Spectrum
	SpecularReflect(ctx context.Context, ray *RayDifferential, si *SurfaceInteraction, scene Scene, sampler Sampler, depth int) Spectrum
	SpecularTransmit(ctx context.Context, ray *RayDifferential, si *SurfaceInteraction, scene Scene, sampler Sampler, depth int) Spectrum
}

func UniformSampleAllLights(it Interaction, scene Scene, sampler Sampler, nLightSamples []int32, handleMedia bool) Spectrum {
	L := NewSpectrum(0)
	for j, light := range scene.Lights() {
		// accumulate contribution of jth light to L
		//light := scene.lights[j]
		nSamples := nLightSamples[j]
		uLightSlice := sampler.Get2DArray(nSamples)
		uScatteringSlice := sampler.Get2DArray(nSamples)
		if len(uLightSlice) == 0 || len(uScatteringSlice) == 0 {
			// use a single sample for illumination from light
			uLight := sampler.Get2D()
			uScattering := sampler.Get2D()
			L.AddAssign(EstimateDirect(it, *uScattering, light, *uLight, scene, sampler, handleMedia, false))
		} else {
			Ld := NewSpectrum(0)
			for k := int32(0); k < nSamples; k++ {
				Ld.AddAssign(EstimateDirect(it, uScatteringSlice[k], light, uLightSlice[k], scene, sampler, handleMedia, false))
			}

			L.AddAssign(Ld.DivScalar(float64(nSamples)))
		}
	}
	return L
}

func UniformSampleOneLight(it Interaction, scene Scene, sampler Sampler, handleMedia bool, lightDistrib *Distribution1D) Spectrum {
	// randomly choose a single light to sample
	nLights := len(scene.Lights())
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

	light := scene.Light(lightNum)
	uLight := sampler.Get2D()
	uScattering := sampler.Get2D()

	spectrum := EstimateDirect(it, *uScattering, light, *uLight, scene, sampler, handleMedia, false)
	spectrum.DivScalar(lightPdf)
	return spectrum
}

func EstimateDirect(it Interaction, uScattering Point2f, light Light, uLight Point2f, scene Scene, sampler Sampler, handleMedia bool, specular bool) Spectrum {
	bsdfFlags := BSDFAll
	if !specular {
		bsdfFlags = BxDFType(BSDFAll &^ BSDFSpecular)
	}

	scatteringPdf := 0.0
	Ld := NewSpectrum(0)
	Li, wi, lightPdf, visibility := light.SampleLi(it, &uLight)

	if lightPdf > 0 && !Li.IsBlack() {
		// compute BSDF for light sampling strategy
		var f Spectrum

		switch isect := it.(type) {
		case *SurfaceInteraction:
			f = isect.BSDF.F(isect.Wo, wi, bsdfFlags)
			wiDotNormal := wi.AbsDot(isect.shading.normal)
			//return NewSpectrum(wiDotNormal)
			f = f.MulScalar(wiDotNormal)
			scatteringPdf = isect.BSDF.Pdf(isect.Wo, wi, bsdfFlags)
		case *MediumInteraction:
			p := isect.phase.P(isect.Wo, wi)
			f = NewSpectrum(p)
			scatteringPdf = p
		default:
			log.Panicf("UnknownInteraction: %+v\n", isect)
		}

		if !f.IsBlack() {
			// compute effect of visibility for light source sample
			if handleMedia {
				Li.MulAssign(visibility.Tr(scene, sampler))
			} else {
				if !visibility.Unoccluded(scene) {
					// Shadow ray blocked
					Li = NewSpectrum(0)
				}
			}

			// add light's contribution to reflected radiance
			if !Li.IsBlack() {
				//fmt.Println("Li not black", IsDeltaLight(light.GetFlags()))
				if IsDeltaLight(light.GetFlags()) {
					Ld.AddAssign(f.Mul(Li).DivScalar(lightPdf))
				} else {
					weight := PowerHeuristic(1.0, lightPdf, 1.0, scatteringPdf)
					Ld.AddAssign(f.Mul(Li).MulScalar(weight).DivScalar(lightPdf))
				}
			}
		}
	}

	// sample SDF with multiple importance sampling
	if !IsDeltaLight(light.GetFlags()) {
		var f Spectrum
		sampledSpecular := false
		switch intr := it.(type) {
		case *SurfaceInteraction:
			// sample scattered Direction for surface interactions
			var sampledType BxDFType
			f, wi, scatteringPdf, sampledType = intr.BSDF.SampleF(intr.Wo, &uScattering, bsdfFlags)
			f.MulScalar(wi.AbsDot(intr.shading.normal))
			sampledSpecular = sampledType&BSDFSpecular != 0
		case *MediumInteraction:
			// sample scattered Direction for Medium interactions
			p := intr.phase.SampleP(intr.Wo, wi, &uScattering)
			f = NewSpectrum(p)
			scatteringPdf = p
		default:
			log.Panic(intr)
		}

		if !f.IsBlack() && scatteringPdf > 0.0 {

			// Account for light contributions along sampled Direction wi
			weight := 1.0
			if !sampledSpecular {
				lightPdf = light.PdfLi(it, wi)
				if lightPdf == 0 {
					return Ld
				}
				weight = PowerHeuristic(1.0, scatteringPdf, 1.0, lightPdf)
			}

			// find intersection and compute transmittance
			ray := it.SpawnRay(wi)
			Tr := NewSpectrum(1.0)
			var foundSurfaceInteraction bool
			lightSI := NewSurfaceInteraction()
			if handleMedia {
				foundSurfaceInteraction = scene.IntersectTr(ray, lightSI, sampler, Tr)
			} else {
				foundSurfaceInteraction = scene.Intersect(ray, lightSI)
			}

			// add light contribution from material sampling
			Li := NewSpectrum(0)
			if foundSurfaceInteraction {
				if lightSI.Primitive.GetAreaLight() == light {
					Li = lightSI.Le(wi.MulScalar(-1))
				}
			} else {
				Li = light.Le(NewRayDifferentialFromRay(ray))
			}
			if !Li.IsBlack() {
				Ld.AddAssign(
					f.Mul(Li).
						Mul(Tr).
						MulScalar(weight / scatteringPdf),
				)
			}
		}
	}

	return Ld
}

type SamplerIntegrator struct {
	camera      Camera
	sampler     Sampler
	pixelBounds *Bounds2i
}

func NewSamplerIntegrator(camera Camera, sampler Sampler, pixelBounds *Bounds2i) *SamplerIntegrator {
	return &SamplerIntegrator{
		camera:      camera,
		sampler:     sampler,
		pixelBounds: pixelBounds,
	}
}

func (s *SamplerIntegrator) GetSampler() Sampler {
	return s.sampler
}

func (s *SamplerIntegrator) GetCamera() Camera {
	return s.camera
}

func (s *SamplerIntegrator) Preprocess(scene Scene, sampler Sampler) {

}

type RenderableTile struct {
	Sampler Sampler
	Bounds  Bounds2i
}

func renderWorker(ctx context.Context, s Integrator, scene Scene, tiles <-chan RenderableTile, progress ProgressReporter) func() error {
	return func() error {
		for rt := range tiles {
			camera := s.GetCamera()

			filmTile := camera.GetFilm().GetFilmTile(&rt.Bounds)

			for pixelY := rt.Bounds.Min.Y; pixelY < rt.Bounds.Max.Y; pixelY++ {
				for pixelX := rt.Bounds.Min.X; pixelX < rt.Bounds.Max.X; pixelX++ {
					pixelPos := &Point2i{pixelX, pixelY}
					rt.Sampler.StartPixel(pixelPos)

					for rt.Sampler.StartNextSample() {
						// initialize CameraSample for current sample
						cameraSample := rt.Sampler.GetCameraSample(pixelPos)

						// generate camera ray for current sample

						rayWeight, rd := camera.GenerateRayDifferential(cameraSample)
						rd.ScaleDifferentials(1.0 / math.Sqrt(float64(rt.Sampler.GetSamplesPerPixel())))

						// evaluate radiance along camera ray
						L := NewSpectrum(0)
						if rayWeight > 0 {
							L = s.Li(ctx, rd, scene, rt.Sampler, 0)
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
						filmTile.AddSample(cameraSample.pFilm, L, rayWeight)

					}

					//px := camera.GetFilm().getPixel(pixelPos)
					//px.value = [3]float64{
					//	L.Index(0),
					//	L.Index(1),
					//	L.Index(2),
					//}

					// merge image tile into film
					//camera.GetFilm().MergeFilmTile(filmTile)

				}
			}

			camera.GetFilm().MergeFilmTile(filmTile)

			progress.Step()
		}

		return nil
	}
}

func Render(ctx context.Context, s Integrator, scene Scene, tileSize int64) error {
	s.Preprocess(scene, s.GetSampler())
	// render image tiles in parallel
	camera := s.GetCamera()

	// compute number of tiles to use for parallel computing
	sampleBounds := camera.GetFilm().GetSampleBounds()
	sampleExtent := sampleBounds.Diagonal()
	nTiles := &Point2i{X: (sampleExtent.X + tileSize - 1) / tileSize, Y: (sampleExtent.Y + tileSize - 1) / tileSize}

	tiles := make(chan RenderableTile)
	progress := NewProgress(nTiles.X * nTiles.Y)
	progress.Start(ctx)

	g, gtx := errgroup.WithContext(ctx)

	for i := 0; i < 64; i++ {
		g.Go(renderWorker(gtx, s, scene, tiles, progress))
	}

	g.Go(func() error {
		defer close(tiles)

		var tileY, tileX int64
		for tileY = 0; tileY < nTiles.Y; tileY++ {
			for tileX = 0; tileX < nTiles.X; tileX++ {
				tile := Point2i{tileX, tileY}
				seed := uint64(tile.Y*nTiles.X + tile.X)

				// compute sample bounds for tile
				x0 := sampleBounds.Min.X + tile.X*tileSize
				x1 := int64(math.Min(float64(x0+tileSize), float64(sampleBounds.Max.X)))
				y0 := sampleBounds.Min.Y + tile.Y*tileSize
				y1 := int64(math.Min(float64(y0+tileSize), float64(sampleBounds.Max.Y)))
				bounds := Bounds2i{&Point2i{x0, y0}, &Point2i{x1, y1}}

				renderableTile := RenderableTile{
					Sampler: s.GetSampler().Clone(seed),
					Bounds:  bounds,
				}

				select {
				case <-gtx.Done():
					return gtx.Err()
				case tiles <- renderableTile:
				}
			}
		}
		return nil
	})

	err := g.Wait()
	if err != nil {
		return errors.Wrap(err, "waiting for render group")
	}

	camera.GetFilm().WriteImage(1.0)

	return nil
}

func SamplerIntegratorSpecularReflect(s Integrator, ctx context.Context, ray *RayDifferential, si *SurfaceInteraction, scene Scene, sampler Sampler, depth int) Spectrum {
	wo := si.Wo
	t := BxDFType(BSDFReflection | BSDFSpecular)
	f, wi, pdf, t := si.BSDF.SampleF(wo, sampler.Get2D(), t)

	// return contribution of specular reflection
	ns := si.shading.normal
	if pdf > 0 && !f.IsBlack() && wi.AbsDot(ns) != 0.0 {
		// compute ray differential rd for specular reflection
		rd := NewRayDifferentialFromRay(si.SpawnRay(wi))
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = si.Point.Add(si.dpdx)
			rd.ryOrigin = si.Point.Add(si.dpdy)

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
		return f.Mul(s.Li(ctx, rd, scene, sampler, depth+1)).MulScalar(wi.AbsDot(ns) / pdf)
	}

	return NewSpectrum(0)
}

func SamplerIntegratorSpecularTransmit(s Integrator, ctx context.Context, ray *RayDifferential, si *SurfaceInteraction, scene Scene, sampler Sampler, depth int) Spectrum {
	wo := si.Wo
	f, wi, pdf, _ := si.BSDF.SampleF(wo, sampler.Get2D(), BxDFType(BSDFTransmission|BSDFSpecular))

	p := si.Point
	ns := si.shading.normal

	L := NewSpectrum(0)
	if pdf > 0 && !f.IsBlack() && wi.AbsDot(ns) != 0.0 {
		rd := NewRayDifferentialFromRay(si.SpawnRay(wi))
		if ray.hasDifferentials {
			rd.hasDifferentials = true
			rd.rxOrigin = p.Add(si.dpdx)
			rd.ryOrigin = p.Add(si.dpdy)

			eta := si.BSDF.Eta
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
		L = f.Mul(s.Li(ctx, rd, scene, sampler, depth+1)).MulScalar(wi.AbsDot(ns) / pdf)
	}

	return L
}
