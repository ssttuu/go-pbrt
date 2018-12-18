package main

import (
	"context"
	"net"

	"google.golang.org/grpc"
)

func serverFunc(ctx context.Context, s *grpc.Server, lis net.Listener) func() error {
	return func() error {
		errch := make(chan error, 1)
		go func() {
			defer close(errch)
			errch <- s.Serve(lis)
		}()

		select {
		case <-ctx.Done():
			s.Stop()
			return ctx.Err()
		case err := <-errch:
			return err
		}
	}
}
