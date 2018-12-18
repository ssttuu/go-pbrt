package main

import (
	"context"
	"fmt"
	"log"
	"net"

	"github.com/ssttuu/go-pbrt/internal/render"
	"github.com/ssttuu/go-pbrt/internal/signal"
	"golang.org/x/sync/errgroup"
	"google.golang.org/grpc"
	"google.golang.org/grpc/reflection"
)

func main() {
	lis, err := net.Listen("tcp", ":3001")
	if err != nil {
		log.Fatalf("failed to listen: %v", err)
	}

	s := grpc.NewServer()
	err = render.RegisterServer(s)
	if err != nil {
		log.Fatalf("failed to register service: %v", err)
	}

	reflection.Register(s)
	fmt.Printf("Listening on %v", lis.Addr())

	g, gtx := errgroup.WithContext(context.Background())
	g.Go(signal.Func(gtx))
	g.Go(serverFunc(gtx, s, lis))

	if err := g.Wait(); err != nil {
		log.Fatal(err)
	}
}
