package main

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

func trace(r *pbrt.Ray, shapes []pbrt.Shaper, depth int) *pbrt.Vector3f {
	//var sphere pbrt.Sphere
	foundSphere := false

	// Find sphere intersection
	for i := 0; i < len(shapes); i++ {
		intersects := shapes[i].IntersectP(r, true)
		if intersects {
			foundSphere = true
		}
	}

	if !foundSphere {
		return &pbrt.Vector3f{0.1, 0.1, 0.1}
	}

	return &pbrt.Vector3f{1, 1, 1}
}

func worker(parentSpan opentracing.Span, id int, camera CameraSettings, shapes []pbrt.Shaper, traceQueue <-chan *TraceData, traceResults chan <- *TraceResult, done chan <- bool) {
	span := parentSpan.Tracer().StartSpan("worker",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer span.Finish()

	//previousSpanContext := span.Context()

	fmt.Println("Worker", id, "started")
	itemProcessing := 0

	var xx, yy float64

	for traceItem := range traceQueue {
		xx = (2.0 * ((float64(traceItem.x) + 0.5) * camera.InverseWidth) - 1.0) * camera.Angle * camera.AspectRatio
		yy = (1.0 - 2.0 * ((float64(traceItem.y) + 0.5) * camera.InverseHeight)) * camera.Angle

		ray := pbrt.NewRay(&pbrt.Point3f{0, 0, 0}, new(pbrt.Vector3f).Set(xx, yy, -1).Normalized(), 0.0)

		pixel := trace(ray, shapes, traceItem.depth)
		traceResults <- &TraceResult{
			x: traceItem.x,
			y: traceItem.y,
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
		Width: width,
		Height: height,
		InverseWidth: 1 / float64(width),
		InverseHeight: 1 / float64(width),
		FOV: fov,
		AspectRatio: float64(width) / float64(height),
		Angle: math.Tan(math.Pi * 0.5 * fov / 180.0),

	}
}

func render(parentSpan opentracing.Span, shapes []pbrt.Shaper) {
	renderSpan := parentSpan.Tracer().StartSpan("render",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer renderSpan.Finish()

	camera := NewCameraSettings(192, 108)
	//camera := NewCameraSettings(1920, 1080)

	imgPng := image.NewNRGBA(image.Rect(0, 0, camera.Width, camera.Height))

	traceQueue := make(chan *TraceData, 2048)
	traceResults := make(chan *TraceResult, 2048)
	nWorkers := 8
	done := make(chan bool, nWorkers)
	for i := 0; i < nWorkers; i++ {
		go worker(renderSpan, i, camera, shapes, traceQueue, traceResults, done)
	}

	go func() {
		populateTraceQueue := parentSpan.Tracer().StartSpan("populate-trace",
			opentracing.ChildOf(renderSpan.Context()),
		)
		for y := 0; y < camera.Height; y++ {
			for x := 0; x < camera.Width; x++ {

				traceQueue <- &TraceData{
					x: x,
					y: y,
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
			Type: "const",
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

	var shapes []pbrt.Shaper

	n := 9

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				x := (float64(i) / float64(n) * 200) - 100
				y := (float64(j) / float64(n) * 200) - 100
				z := (float64(k) / float64(n) * 200) - 100
				radius := 1.0
				if x > 0 {
					radius = 2.0
				}
				if z < 0 {
					radius *= 2
				}

				xform := pbrt.Translate(&pbrt.Vector3f{x, y, z})

				shapes = append(shapes, pbrt.NewSphereShape(xform, xform.Inverse(), false, radius))
			}
		}
	}

	span := opentracing.StartSpan("render")
	defer span.Finish()

	render(span, shapes)
}
