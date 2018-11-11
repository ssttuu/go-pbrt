package render

import (
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
	"math"
	"image"
	"image/color"
	"os"
	"image/png"
	"log"
	"fmt"

	"github.com/uber/jaeger-client-go/config"
	"github.com/uber/jaeger-lib/metrics/prometheus"
	"github.com/opentracing/opentracing-go"
	"io"
	"flag"
	"runtime/pprof"
)

type TraceData struct {
	x, y  int
	depth int
}

type TraceResult struct {
	x, y  int
	value pbrt.Vector3f
}

func trace(r *pbrt.Ray, scene *pbrt.Scene, depth int) *pbrt.Vector3f {
	foundSphere := scene.Aggregate.IntersectP(r)

	if !foundSphere {
		return &pbrt.Vector3f{0.1, 0.1, 0.1}
	}

	return &pbrt.Vector3f{1, 1, 1}
}

func worker(parentSpan opentracing.Span, id int, camera CameraSettings, scene *pbrt.Scene, traceQueue <-chan *TraceData, traceResults chan<- *TraceResult, done chan<- bool) {
	span := parentSpan.Tracer().StartSpan("worker",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer span.Finish()

	//previousSpanContext := span.Context()

	fmt.Println("Worker", id, "started")
	itemProcessing := 0

	var xx, yy float64

	for traceItem := range traceQueue {
		xx = (2.0*((float64(traceItem.x)+0.5)*camera.InverseWidth) - 1.0) * camera.Angle * camera.AspectRatio
		yy = (1.0 - 2.0*((float64(traceItem.y)+0.5)*camera.InverseHeight)) * camera.Angle

		ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, new(pbrt.Vector3f).Set(xx, yy, -1).Normalized(), 0.0)

		pixel := trace(ray, scene, traceItem.depth)
		traceResults <- &TraceResult{
			x:     traceItem.x,
			y:     traceItem.y,
			value: *pixel,
		}
		itemProcessing += 1
	}
	done <- true
}

type CameraSettings struct {
	Width         int
	Height        int
	InverseWidth  float64
	InverseHeight float64
	FOV           float64
	AspectRatio   float64
	Angle         float64
}

func NewCameraSettings(width, height int) CameraSettings {
	fov := 30.0
	return CameraSettings{
		Width:         width,
		Height:        height,
		InverseWidth:  1 / float64(width),
		InverseHeight: 1 / float64(width),
		FOV:           fov,
		AspectRatio:   float64(width) / float64(height),
		Angle:         math.Tan(math.Pi * 0.5 * fov / 180.0),
	}
}

func render(parentSpan opentracing.Span, scene *pbrt.Scene) {
	renderSpan := parentSpan.Tracer().StartSpan("render",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer renderSpan.Finish()

	camera := NewCameraSettings(192*5, 108*5)
	//camera := NewCameraSettings(1920, 1080)

	imgPng := image.NewNRGBA(image.Rect(0, 0, camera.Width, camera.Height))

	traceQueue := make(chan *TraceData, 2048)
	traceResults := make(chan *TraceResult, 2048)
	nWorkers := 8
	done := make(chan bool, nWorkers)
	for i := 0; i < nWorkers; i++ {
		go worker(renderSpan, i, camera, scene, traceQueue, traceResults, done)
	}

	go func() {
		populateTraceQueue := parentSpan.Tracer().StartSpan("populate-trace",
			opentracing.ChildOf(renderSpan.Context()),
		)
		for y := 0; y < camera.Height; y++ {
			for x := 0; x < camera.Width; x++ {

				traceQueue <- &TraceData{
					x:     x,
					y:     y,
					depth: 0,
				}
			}
		}
		close(traceQueue)
		populateTraceQueue.Finish()
	}()

	writeImageTrace := parentSpan.Tracer().StartSpan("write-image",
		opentracing.ChildOf(renderSpan.Context()),
	)
	nWorkersCompleted := 0
	for {
		if nWorkersCompleted == nWorkers {
			break
		}

		select {
		case traceResult := <-traceResults:
			imgPng.Set(traceResult.x, traceResult.y, color.NRGBA{
				R: uint8(traceResult.value.X * 255),
				G: uint8(traceResult.value.Y * 255),
				B: uint8(traceResult.value.Z * 255),
				A: 255,
			})
		case <-done:
			nWorkersCompleted += 1
		}
	}

	close(traceResults)

	f, err := os.Create("image.png")
	if err != nil {
		log.Fatal(err)
	}

	if err := png.Encode(f, imgPng); err != nil {
		f.Close()
		log.Fatal(err)
	}

	if err := f.Close(); err != nil {
		log.Fatal(err)
	}

	writeImageTrace.Finish()

}

func getTracer() (opentracing.Tracer, io.Closer) {
	metricsFactory := prometheus.New()
	cfg := config.Configuration{
		Sampler: &config.SamplerConfig{
			Type:  "const",
			Param: 1,
		},
		Reporter: &config.ReporterConfig{
			LocalAgentHostPort: "jaeger:6831",
		},
	}
	tracer, closer, err := cfg.New(
		"pbrt",
		config.Metrics(metricsFactory),
	)
	if err != nil {
		log.Fatal("Failed to start tracer", err)
	}

	return tracer, closer

}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	//tracer, closer := getTracer()
	//defer closer.Close()

	// START

	var primitives []pbrt.Primitiver

	n := 8

	kd := pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(0.5))
	sigma := pbrt.NewConstantFloatTexture(0.0)
	material := pbrt.NewMatteMaterial(kd, sigma, nil)

	for i := 0; i <= n; i++ {
		for j := 0; j <= n; j++ {
			for k := 0; k <= n; k++ {
				x := (float64(i) / float64(n) * 200) - 100
				y := (float64(j) / float64(n) * 200) - 100
				z := (float64(k) / float64(n) * 200) - 100
				radius := 4.0
				//if x > 0 {
				//	radius = 1.0
				//}
				//if z < 0 {
				//	radius *= 2
				//}

				xform := pbrt.Translate(&pbrt.Vector3f{x, y, z})
				sphere := pbrt.NewSphereShape(
					fmt.Sprintf("Sphere: %v, %v, %v - MatteMaterial", x, y, z),
					xform, xform.Inverse(), false, radius)
				geoPrim := pbrt.NewGeometricPrimitive(sphere, material)

				primitives = append(primitives, geoPrim)
			}
		}
	}

	span := opentracing.StartSpan("render")
	defer span.Finish()

	agg := pbrt.NewSimpleAggregate(primitives)

	//areaLightXform := pbrt.Translate(&pbrt.Vector3f{0, 0, 0})
	//areaLightShape := pbrt.NewSphereShape(areaLightXform, areaLightXform.Inverse(), false, 10.0)
	//areaLight := pbrt.NewDiffuseAreaLight(areaLightXform, nil, pbrt.NewSpectrum(100000.0), 10, areaLightShape, false)

	pointLight := pbrt.NewPointLight(pbrt.Translate(&pbrt.Vector3f{-20, 0, 0}), nil, pbrt.NewSpectrum(pbrt.Pi * 100000))

	lights := []pbrt.Lighter{
		pointLight,
		//areaLight,
	}

	scene := &pbrt.Scene{
		Aggregate:  agg,
		Lights: lights,
		WorldBound: agg.WorldBound(),
	}

	camXform := pbrt.Translate(&pbrt.Vector3f{10, 10, 0}).Mul(pbrt.RotateY(-45))
	camAnimXform := pbrt.NewAnimatedTransform(camXform, camXform, 0, 1)

	shutterOpen, shutterClose := 0.0, 1.0
	div := int64(1)
	resolution := &pbrt.Point2i{1920 / div, 1080 / div}
	cropBounds := &pbrt.Bounds2f{Min: &pbrt.Point2f{0, 0}, Max: &pbrt.Point2f{1, 1}}
	boxFilter := pbrt.NewBoxFilter(&pbrt.Point2f{1, 1})
	film := pbrt.NewFilm("image-new.png", resolution, cropBounds, boxFilter, 100.0, 1.0, 1.0)
	camera := pbrt.NewPerspectiveCamera(camAnimXform, cropBounds, shutterOpen, shutterClose, 0, 20, 90.0, film, nil)

	pixelBounds := &pbrt.Bounds2i{&pbrt.Point2i{0, 0}, resolution}
	sampler := pbrt.NewRandomSampler(10, 4)
	integrator := pbrt.NewDirectLightingIntegrator(pbrt.UniformSampleAll, 10, camera, sampler, pixelBounds)
	pbrt.Render(integrator, scene)

	//render(span, scene)
}
