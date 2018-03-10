package pbrt

type Cameraer interface {

}

type Camera struct {
	cameraToWorld             *AnimatedTransform
	shutterOpen, shutterClose float64
	Film                      *Film
	medium                    Mediumer
}

func NewCamera(cameraToWorld *AnimatedTransform, shutterOpen, shutterClose float64, film *Film, medium Mediumer) *Camera {
	return &Camera{
		cameraToWorld: cameraToWorld,
		shutterOpen:   shutterOpen,
		shutterClose:  shutterClose,
		Film:          film,
		medium:        medium,
	}
}

func (c *Camera) GenerateRay(sample *CameraSample) (float64, *Ray) {
	return 0, nil
}

func (c *Camera) GenerateRayDifferential(sample *CameraSample) (float64, *RayDifferential) {
	wt, r := c.GenerateRay(sample)
	rd := &RayDifferential{
		Ray: r,
	}

	// find camera ray after shifting one pixel in the x direction
	sshift := sample
	sshift.pFilm.X++

	wtx, rx := c.GenerateRay(sshift)
	if wtx == 0.0 {
		return 0.0, rd
	}
	rd.rxOrigin = rx.origin
	rd.rxDirection = rx.direction

	// find camera ray after shifting one pixel in the y direction
	sshift.pFilm.X--
	sshift.pFilm.Y++

	wty, ry := c.GenerateRay(sshift)
	if wty == 0.0 {
		return 0.0, rd
	}

	rd.ryOrigin = ry.origin
	rd.ryDirection = ry.direction
	rd.hasDifferentials = true
	return wt, rd
}

func (c *Camera) We(r *Ray, pRaster *Point2f) Spectrum {
	return nil
}

func (c *Camera) PdfWe(r *Ray, pdfPos, pdfDir float64) {

}

func (c *Camera) SampleWi(ref Interaction, u *Point2f, wi *Vector3f, pdf float64, pRaser *Point2f, vis *VisibilityTester) Spectrum {
	return nil
}

type CameraSample struct {
	pFilm *Point2f
	pLens *Point2f
	time  float64
}

type ProjectiveCamera struct {
	*Camera

	CameraToScreen, RasterToCamera *Transform
	ScreenToRaster, RasterToScreen *Transform
	lensRadius, focalDistance      float64
}

func NewProjectiveCamera(cameraToWorld *AnimatedTransform, cameraToScreen *Transform, screenWindow *Bounds2f, shutterOpen, shutterClose float64, lensr, focald float64, film *Film, medium Mediumer) *ProjectiveCamera {

	screenToRaster := Scale(float64(film.FullResolution.X), float64(film.FullResolution.Y), 0)
	screenToRaster = screenToRaster.Mul(Scale(1.0 / (screenWindow.max.X - screenWindow.min.X), 1.0 / (screenWindow.min.Y - screenWindow.max.Y), 1))
	screenToRaster = screenToRaster.Mul(Translate(&Vector3f{-screenWindow.min.X, -screenWindow.max.Y, 0}))

	rasterToScreen := screenToRaster.Inverse()
	rasterToCamera := cameraToScreen.Inverse().Mul(rasterToScreen)


	return &ProjectiveCamera{
		Camera: NewCamera(cameraToWorld, shutterOpen, shutterOpen, film, medium),
		CameraToScreen: cameraToScreen,
		RasterToCamera: rasterToCamera,
		ScreenToRaster: screenToRaster,
		RasterToScreen: rasterToScreen,
		lensRadius: lensr,
		focalDistance: focald,
	}
}
