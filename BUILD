cc_binary(
	name = "1d-images-solution-likelihood",
	srcs = ["1d-images-solution-likelihood.cpp"],
	deps = [":1d-advection-diffusion-images"],
)

cc_binary(
	name = "2d-images-solution-test",
	srcs = ["2d-images-solution-test.cpp"],
	deps = [":2d-advection-diffusion-images"],
)

cc_library(
	name = "1d-advection-diffusion-images",
	srcs = ["1DAdvectionDiffusionSolverImages.cpp"],
	hdrs = ["1DAdvectionDiffusionSolverImages.hpp"],
	visibility = ["//visibility:public"],
	linkopts = ["-lm", "-lgsl", "-lgslcblas"]
)

cc_library(
	name = "2d-advection-diffusion-images",
	srcs = ["2DAdvectionDiffusionSolverImages.cpp"],
	hdrs = ["2DAdvectionDiffusionSolverImages.hpp"],
	copts = ["-Isrc/images-expansion"],
	linkopts = ["-lm", "-lgsl", "-lgslcblas"],
	deps = [":1d-advection-diffusion-images"],
	visibility = ["//visibility:public"],
)

cc_library(
	name = "cmath",
	linkopts = ["-lm"]
)

cc_library(
	name = "gsl",
	linkopts = ["-lgsl", "-lgslcblas", "-lm"],
)

cc_library(
	name = "gslcblas",
	linkopts = ["-lgslcblas"],
)