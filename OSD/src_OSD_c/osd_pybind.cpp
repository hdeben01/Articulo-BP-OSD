//DESCOMENTAR CÓDIGO DE ABAJO Y EJECUTAR(Linux) :
//        Se debe tener que cmake --version de ~3.28. Para ello, desde GaBP_cod_cpp ejecutar:
//
//        export PATH=$HOME/.local/bin:$PATH
//
//        De esta forma se coge la versión local del cmake, que ya está descargada.
//
//        Cuando se tenga esta versión del cmake, compilar con:
//
//        cmake -Dpybind11_DIR=$(python3 -m pybind11 --cmakedir) ..
//
//        make -j$(nproc)
//        *** IMPORTANTE ***
//        Cuando se tenga el .so, se debe de ejecutar en un entorno Python3.13


/*
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <memory>
#include "OSD_python.h"

namespace py = pybind11;

PYBIND11_MODULE(OSD_decoder, m) {
    m.doc() = "OSD decoder";

    py::class_<osd_decoder>(m, "OSD_decoder")
        .def(py::init([](py::array_t<int, py::array::c_style | py::array::forcecast> H_cod_in) {
            auto buf = H_cod_in.request();
            int* H_ptr = static_cast<int*>(buf.ptr);

            // Instanciación directa en Heap para no agotar el Stack
            auto decoder = std::make_unique<osd_decoder>();

            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M; j++) {
                    decoder->H[i][j] = H_ptr[i * M + j];
                }
            }

            return decoder.release();
        }), py::arg("H_cod"))

        .def("decode", [](osd_decoder& self, 
                          py::array_t<int, py::array::c_style | py::array::forcecast> synd_in, 
                          py::array_t<double, py::array::c_style | py::array::forcecast> prob_ini_in) {

            auto synd_buf = synd_in.request();
            int* synd_ptr = static_cast<int*>(synd_buf.ptr);

            auto prob_buf = prob_ini_in.request();
            double* prob_ptr = static_cast<double*>(prob_buf.ptr);

            int synd_c[N];
            double prob_ini_c[M];
            int sol_c[M] = { 0 };

            std::copy(synd_ptr, synd_ptr + N, synd_c);
            std::copy(prob_ptr, prob_ptr + M, prob_ini_c);

            self.decode(synd_c, prob_ini_c, sol_c);

            py::array_t<int> result(M);
            auto res_buf = result.request();
            int* res_ptr = static_cast<int*>(res_buf.ptr);

            std::copy(sol_c, sol_c + M, res_ptr);

            return result;
        }, py::arg("synd"), py::arg("prob_ini"));
}*/