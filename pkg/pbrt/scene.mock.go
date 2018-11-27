// Code generated by MockGen. DO NOT EDIT.
// Source: scene.go

// Package pbrt is a generated GoMock package.
package pbrt

import (
	gomock "github.com/golang/mock/gomock"
	reflect "reflect"
)

// MockScene is a mock of Scene interface
type MockScene struct {
	ctrl     *gomock.Controller
	recorder *MockSceneMockRecorder
}

// MockSceneMockRecorder is the mock recorder for MockScene
type MockSceneMockRecorder struct {
	mock *MockScene
}

// NewMockScene creates a new mock instance
func NewMockScene(ctrl *gomock.Controller) *MockScene {
	mock := &MockScene{ctrl: ctrl}
	mock.recorder = &MockSceneMockRecorder{mock}
	return mock
}

// EXPECT returns an object that allows the caller to indicate expected use
func (m *MockScene) EXPECT() *MockSceneMockRecorder {
	return m.recorder
}

// Aggregate mocks base method
func (m *MockScene) Aggregate() Aggregate {
	ret := m.ctrl.Call(m, "Aggregate")
	ret0, _ := ret[0].(Aggregate)
	return ret0
}

// Aggregate indicates an expected call of Aggregate
func (mr *MockSceneMockRecorder) Aggregate() *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "Aggregate", reflect.TypeOf((*MockScene)(nil).Aggregate))
}

// Intersect mocks base method
func (m *MockScene) Intersect(r *Ray, si *SurfaceInteraction) bool {
	ret := m.ctrl.Call(m, "Intersect", r, si)
	ret0, _ := ret[0].(bool)
	return ret0
}

// Intersect indicates an expected call of Intersect
func (mr *MockSceneMockRecorder) Intersect(r, si interface{}) *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "Intersect", reflect.TypeOf((*MockScene)(nil).Intersect), r, si)
}

// IntersectP mocks base method
func (m *MockScene) IntersectP(r *Ray) bool {
	ret := m.ctrl.Call(m, "IntersectP", r)
	ret0, _ := ret[0].(bool)
	return ret0
}

// IntersectP indicates an expected call of IntersectP
func (mr *MockSceneMockRecorder) IntersectP(r interface{}) *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "IntersectP", reflect.TypeOf((*MockScene)(nil).IntersectP), r)
}

// IntersectTr mocks base method
func (m *MockScene) IntersectTr(r *Ray, si *SurfaceInteraction, sampler Sampler, transmittance Spectrum) bool {
	ret := m.ctrl.Call(m, "IntersectTr", r, si, sampler, transmittance)
	ret0, _ := ret[0].(bool)
	return ret0
}

// IntersectTr indicates an expected call of IntersectTr
func (mr *MockSceneMockRecorder) IntersectTr(r, si, sampler, transmittance interface{}) *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "IntersectTr", reflect.TypeOf((*MockScene)(nil).IntersectTr), r, si, sampler, transmittance)
}

// Light mocks base method
func (m *MockScene) Light(index int) Light {
	ret := m.ctrl.Call(m, "Light", index)
	ret0, _ := ret[0].(Light)
	return ret0
}

// Light indicates an expected call of Light
func (mr *MockSceneMockRecorder) Light(index interface{}) *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "Light", reflect.TypeOf((*MockScene)(nil).Light), index)
}

// Lights mocks base method
func (m *MockScene) Lights() []Light {
	ret := m.ctrl.Call(m, "Lights")
	ret0, _ := ret[0].([]Light)
	return ret0
}

// Lights indicates an expected call of Lights
func (mr *MockSceneMockRecorder) Lights() *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "Lights", reflect.TypeOf((*MockScene)(nil).Lights))
}

// InfiniteLights mocks base method
func (m *MockScene) InfiniteLights() []Light {
	ret := m.ctrl.Call(m, "InfiniteLights")
	ret0, _ := ret[0].([]Light)
	return ret0
}

// InfiniteLights indicates an expected call of InfiniteLights
func (mr *MockSceneMockRecorder) InfiniteLights() *gomock.Call {
	return mr.mock.ctrl.RecordCallWithMethodType(mr.mock, "InfiniteLights", reflect.TypeOf((*MockScene)(nil).InfiniteLights))
}
