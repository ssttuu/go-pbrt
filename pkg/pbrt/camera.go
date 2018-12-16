//go:generate mockgen -source=camera.go -destination=camera.mock.go -package=pbrt

package pbrt

import "github.com/ssttuu/go-pbrt/pkg/math"

type Camera interface {
	GenerateRay(sample *CameraSample) (float64, *Ray)
	GenerateRayDifferential(sample *CameraSample) (float64, *RayDifferential)
	We(r *Ray) (Spectrum, *Point2f)
	PdfWe(r *Ray) (pdfPos, pdfDir float64)
	SampleWi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, pRaster *Point2f, vis *VisibilityTester)

	GetFilm() *Film
}

type camera struct {
	cameraToWorld             *AnimatedTransform
	shutterOpen, shutterClose float64
	Film                      *Film
	medium                    Medium
}

func NewCamera(cameraToWorld *AnimatedTransform, shutterOpen, shutterClose float64, film *Film, medium Medium) *camera {
	return &camera{
		cameraToWorld: cameraToWorld,
		shutterOpen:   shutterOpen,
		shutterClose:  shutterClose,
		Film:          film,
		medium:        medium,
	}
}

func (c *camera) GetFilm() *Film {
	return c.Film
}

//func (c *camera) GenerateRay(sample *CameraSample) (float64, *Ray) {
//	return 0, nil
//}
//
//func (c *camera) GenerateRayDifferential(sample *CameraSample) (float64, *RayDifferential) {
//	wt, r := c.GenerateRay(sample)
//	rd := NewRayDifferentialFromRay(r)
//
//	// find camera ray after shifting one pixel in the x Direction
//	sshift := sample
//	sshift.pFilm.X++
//
//	wtx, rx := c.GenerateRay(sshift)
//	if wtx == 0.0 {
//		return 0.0, rd
//	}
//	rd.rxOrigin = rx.Origin
//	rd.rxDirection = rx.Direction
//
//	// find camera ray after shifting one pixel in the y Direction
//	sshift.pFilm.X--
//	sshift.pFilm.Y++
//
//	wty, ry := c.GenerateRay(sshift)
//	if wty == 0.0 {
//		return 0.0, rd
//	}
//
//	rd.ryOrigin = ry.Origin
//	rd.ryDirection = ry.Direction
//	rd.hasDifferentials = true
//	return wt, rd
//}

//func (c *camera) We(r *Ray, pRaster *Point2f) Spectrum {
//	return nil
//}
//
//func (c *camera) PdfWe(r *Ray, pdfPos, pdfDir float64) {
//
//}
//
//func (c *camera) SampleWi(ref Interaction, u *Point2f, wi *Vector3f, pdf float64, pRaser *Point2f, vis *VisibilityTester) Spectrum {
//	return nil
//}

type CameraSample struct {
	pFilm *Point2f
	pLens *Point2f
	time  float64
}

type ProjectiveCamera struct {
	*camera

	CameraToScreen, RasterToCamera *Transform
	ScreenToRaster, RasterToScreen *Transform
	lensRadius, focalDistance      float64
}

func NewProjectiveCamera(cameraToWorld *AnimatedTransform, cameraToScreen *Transform, screenWindow *Bounds2f, shutterOpen, shutterClose float64, lensRadius, focalDistance float64, film *Film, medium Medium) *ProjectiveCamera {

	screenToRaster := Scale(float64(film.FullResolution.X), float64(film.FullResolution.Y), 1.0)
	screenToRaster = screenToRaster.Mul(Scale(1.0/(screenWindow.Max.X-screenWindow.Min.X), 1.0/(screenWindow.Min.Y-screenWindow.Max.Y), 1.0))
	screenToRaster = screenToRaster.Mul(Translate(&Vector3f{-screenWindow.Min.X, -screenWindow.Max.Y, 0}))

	rasterToScreen := screenToRaster.Inverse()
	rasterToCamera := cameraToScreen.Inverse().Mul(rasterToScreen)

	return &ProjectiveCamera{
		camera:         NewCamera(cameraToWorld, shutterOpen, shutterOpen, film, medium),
		CameraToScreen: cameraToScreen,
		RasterToCamera: rasterToCamera,
		ScreenToRaster: screenToRaster,
		RasterToScreen: rasterToScreen,
		lensRadius:     lensRadius,
		focalDistance:  focalDistance,
	}
}

// TODO: move out of pbrt core

type PerspectiveCamera struct {
	*ProjectiveCamera

	dxCamera, dyCamera *Vector3f
	A                  float64
}

func NewPerspectiveCamera(cameraToWorld *AnimatedTransform, screenWindow *Bounds2f, shutterOpen, shutterClose float64, lensRadius, focalDistance, fov float64, film *Film, medium Medium) *PerspectiveCamera {
	projCamera := NewProjectiveCamera(cameraToWorld, Perspective(fov, 1e-2, 1000.0), screenWindow, shutterOpen, shutterClose, lensRadius, focalDistance, film, medium)

	// compute differential changes in Origin for perspective camera rays
	dxCamera := func() *Point3f {
		dxOne, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{1, 0, 0}, NoError)
		dxZero, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{0, 0, 0}, NoError)
		return dxOne.Sub(dxZero)
	}()

	dyCamera := func() *Point3f {
		dyOne, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{0, 1, 0}, NoError)
		dyZero, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{0, 0, 0}, NoError)
		return dyOne.Sub(dyZero)
	}()

	// compute image plane bounds at z=1 for PerspectiveCamera
	res := film.FullResolution
	pMin, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{0, 0, 0}, NoError)
	pMax, _ := projCamera.RasterToCamera.TransformPoint(&Point3f{float64(res.X), float64(res.Y), 0}, NoError)
	pMin = pMin.DivScalar(pMin.Z)
	pMax = pMax.DivScalar(pMax.Z)
	A := math.Abs((pMax.X - pMin.X) * (pMax.Y - pMin.Y))

	return &PerspectiveCamera{
		ProjectiveCamera: projCamera,
		dxCamera:         dxCamera,
		dyCamera:         dyCamera,
		A:                A,
	}
}

func (c *PerspectiveCamera) GenerateRay(sample *CameraSample) (float64, *Ray) {
	// Compute raster and camera sample positions
	pFilm := &Point3f{sample.pFilm.X, sample.pFilm.Y, 0}
	pCamera, _ := c.RasterToCamera.TransformPoint(pFilm, NoError)
	ray := NewRay(&Point3f{0, 0, 0}, pCamera.Normalized(), 0.0)

	// modify ray for depth of field
	if c.lensRadius > 0 {
		// sample Point on lens
		pLens := ConcentricSampleDisk(sample.pLens).MulScalar(c.lensRadius)

		// compute Point on plane of focus
		ft := c.focalDistance / ray.Direction.Z
		pFocus := ray.Direction.MulScalar(ft).Add(ray.Origin)

		// update ray for effect of lens
		ray.Origin = &Point3f{pLens.X, pLens.Y, 0}
		ray.Direction = pFocus.Sub(ray.Origin).Normalized()
	}

	ray.Time = math.Lerp(sample.time, c.shutterOpen, c.shutterClose)
	ray.Medium = c.medium
	return 1.0, c.cameraToWorld.TransformRay(ray)
}

func (c *PerspectiveCamera) GenerateRayDifferential(sample *CameraSample) (float64, *RayDifferential) {
	// compute raster and camera sample position
	pFilm := &Point3f{sample.pFilm.X, sample.pFilm.Y, 0}
	pCamera, _ := c.RasterToCamera.TransformPoint(pFilm, NoError)
	dir := pCamera.Normalized()
	ray := NewRayDifferentialFromRay(NewRay(&Point3f{0, 0, 0}, dir, 0))

	// modify ray for depth of field
	if c.lensRadius > 0 {
		// sample points on lens
		pLens := ConcentricSampleDisk(sample.pLens).MulScalar(c.lensRadius)

		// compute Point on plane of focus
		ft := c.focalDistance / ray.Direction.Z
		pFocus := ray.Direction.MulScalar(ft).Add(ray.Origin)

		// update ray for effect of lens
		ray.Origin = &Point3f{pLens.X, pLens.Y, 0}
		ray.Direction = pFocus.Sub(ray.Origin).Normalized()
	}

	// compute offset rays for PerspectiveCamera ray differentials
	if c.lensRadius > 0 {
		// compute PerspectiveCamera ray differentials accounting for lens

		// sample Point on lens
		pLens := ConcentricSampleDisk(sample.pLens).MulScalar(c.lensRadius)
		dx := pCamera.Add(c.dxCamera).Normalized()
		ft := c.focalDistance / dx.Z
		pFocus := ray.Direction.Mul(dx.MulScalar(ft)).Add(ray.Origin)
		ray.rxOrigin = &Point3f{pLens.X, pLens.Y, 0}
		ray.rxDirection = pFocus.Sub(ray.rxOrigin).Normalized()

		dy := pCamera.Add(c.dyCamera).Normalized()
		ft = c.focalDistance / dy.Z
		pFocus = new(Point3f).Add(dy.MulScalar(ft))
		ray.ryOrigin = &Point3f{pLens.X, pLens.Y, 0}
		ray.ryDirection = pFocus.Sub(ray.ryOrigin).Normalized()
	} else {
		ray.rxOrigin = ray.Origin
		ray.ryOrigin = ray.Origin
		ray.rxDirection = pCamera.Add(c.dxCamera).Normalized()
		ray.ryDirection = pCamera.Add(c.dyCamera).Normalized()
	}
	ray.Time = math.Lerp(sample.time, c.shutterOpen, c.shutterClose)
	ray.Medium = c.medium
	ray.Ray = c.cameraToWorld.TransformRay(ray.Ray)
	ray.hasDifferentials = true

	return 1, ray
}

func (c *PerspectiveCamera) We(r *Ray) (Spectrum, *Point2f) {
	// interpolate camera Matrix and check if w is forward-facing
	c2w := c.cameraToWorld.Interpolate(r.Time)
	cosTheta := r.Direction.Dot(c2w.TransformVector(&Vector3f{0, 0, 1}))
	if cosTheta <= 0 {
		return NewSpectrum(0), nil
	}

	// map ray (p, w) onto the raster grid
	v := 1.0
	if c.lensRadius > 0 {
		v = c.focalDistance
	}
	pFocus := r.Direction.MulScalar(v).Add(r.Origin)
	pRaster := c.RasterToCamera.Inverse().TransformVector(c2w.Inverse().TransformVector(pFocus))

	// return zero imporance for out of bounds points
	sampleBounds := c.Film.GetSampleBounds()
	if pRaster.X < float64(sampleBounds.Min.X) || pRaster.X >= float64(sampleBounds.Max.X) || pRaster.Y < float64(sampleBounds.Min.Y) || pRaster.Y >= float64(sampleBounds.Max.Y) {
		return NewSpectrum(0), nil
	}

	// compute lens area of perspective camera
	lensArea := 1.0
	if c.lensRadius > 0 {
		lensArea = math.Pi * c.lensRadius * c.lensRadius
	}

	// return imporance for Point on image plane
	cos2theta := cosTheta * cosTheta
	return NewSpectrum(1.0 / (c.A * lensArea * cos2theta * cos2theta)), &Point2f{pRaster.X, pRaster.Y}

}

func (c *PerspectiveCamera) PdfWe(r *Ray) (pdfPos, pdfDir float64) {
	c2w := c.cameraToWorld.Interpolate(r.Time)
	cosTheta := r.Direction.Dot(c2w.TransformVector(&Vector3f{0, 0, 0}))
	if cosTheta <= 0 {
		return 0, 0
	}

	// compute lens area of perspective camera
	lensArea := 1.0
	if c.lensRadius > 0 {
		lensArea = math.Pi * c.lensRadius * c.lensRadius
	}

	return 1.0 / lensArea, 1.0 / (c.A * cosTheta * cosTheta * cosTheta)
}

func (c *PerspectiveCamera) SampleWi(ref Interaction, u *Point2f) (s Spectrum, wi *Vector3f, pdf float64, pRaster *Point2f, vis *VisibilityTester) {
	// Uniformly sample a lens interaction
	pLens := ConcentricSampleDisk(u).MulScalar(c.lensRadius)
	pLensWorld := c.cameraToWorld.TransformPointAtTime(&Point3f{pLens.X, pLens.Y, 0}, ref.GetTime())
	lensIntr := &interaction{
		Point:  pLensWorld,
		Normal: c.cameraToWorld.TransformVectorAtTime(&Vector3f{0, 0, 1}, ref.GetTime()),
		time:   ref.GetTime(),
		//mediumAccessor: c.Medium,
	}

	// populate arguments and compute the importance value
	vis = &VisibilityTester{ref, lensIntr}
	wi = lensIntr.Point.Sub(ref.GetPoint())
	dist := wi.Length()
	wi = wi.DivScalar(dist)

	// compute PDF for importance arriving at ref

	// compute lens area of perspective camera
	lensArea := 1.0
	if c.lensRadius > 0 {
		lensArea = math.Pi * c.lensRadius * c.lensRadius
	}

	pdf = (dist * dist) / (lensIntr.Normal.AbsDot(wi) * lensArea)

	spectrum, pRaster := c.We(lensIntr.SpawnRay(wi.MulScalar(-1.0)))

	return spectrum, wi, pdf, pRaster, vis
}
