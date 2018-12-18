
pkg/proto/render:
	mkdir -p pkg/proto/render

pkg/proto/render/service.pb.go: pkg/proto/render
	protoc -I=./proto --go_out=plugins=grpc:./pkg/proto ./proto/render/service.proto

