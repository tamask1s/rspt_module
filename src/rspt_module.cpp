#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Python 2D lista inplace négyzetre emelés
void inplace_square_2dlist(py::list list2d) {
    for (auto row : list2d) {
        for (size_t i = 0; i < py::len(row); ++i) {
            row[i] = py::float_(py::float_(row[i]) * py::float_(row[i]));
        }
    }
}

// Python 2D lista új listába másolva négyzetre emelés
py::list square_2dlist(py::list list2d) {
    py::list result;
    for (auto row : list2d) {
        py::list new_row;
        for (auto value : row) {
            new_row.append(py::float_(py::float_(value) * py::float_(value)));
        }
        result.append(new_row);
    }
    return result;
}

// NumPy tömb inplace négyzetre emelés
void inplace_square_numpy(py::array_t<double> array) {
    auto buf = array.mutable_unchecked<2>();
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        for (ssize_t j = 0; j < buf.shape(1); ++j)
            buf(i, j) *= buf(i, j);
}

// NumPy tömb új példányban négyzetre emelés
py::array_t<double> square_numpy(py::array_t<double> array) {
    auto buf = array.unchecked<2>();
    auto result = py::array_t<double>(array.shape());
    auto res_buf = result.mutable_unchecked<2>();

    for (ssize_t i = 0; i < buf.shape(0); ++i)
        for (ssize_t j = 0; j < buf.shape(1); ++j)
            res_buf(i, j) = buf(i, j) * buf(i, j);

    return result;
}

// Modul definiálás
PYBIND11_MODULE(rspt_module, m) {
    m.doc() = "2D tömbök négyzetre emelése C++ és Python között";

    m.def("inplace_square_2dlist", &inplace_square_2dlist, "2D Python lista inplace négyzetre emelése");
    m.def("square_2dlist", &square_2dlist, "2D Python lista négyzetre emelése új példányban");

    m.def("inplace_square_numpy", &inplace_square_numpy, "NumPy tömb inplace négyzetre emelése");
    m.def("square_numpy", &square_numpy, "NumPy tömb négyzetre emelése új példányban");
}
