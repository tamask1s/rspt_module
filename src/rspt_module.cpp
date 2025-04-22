#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

// Python list 2D inplace square
void square_inplace_list(py::list arr) {
    for (auto row : arr) {
        for (size_t i = 0; i < py::len(row); ++i) {
            row[i] = py::float_(row[i]).cast<double>() * py::float_(row[i]).cast<double>();
        }
    }
}

// Numpy 2D inplace square
void square_inplace_numpy(py::array_t<double> arr) {
    auto buf = arr.mutable_unchecked<2>();
    for (ssize_t i = 0; i < buf.shape(0); i++)
        for (ssize_t j = 0; j < buf.shape(1); j++)
            buf(i, j) *= buf(i, j);
}

// Python list 2D, returns new list
py::list square_return_list(py::list arr) {
    py::list result;
    for (auto row : arr) {
        py::list new_row;
        for (auto val : row) {
            auto v = py::float_(val).cast<double>();
            new_row.append(v * v);
        }
        result.append(new_row);
    }
    return result;
}

// Numpy 2D, returns new numpy array
py::array_t<double> square_return_numpy(py::array_t<double> arr) {
    auto buf = arr.unchecked<2>();
    py::array_t<double> result({buf.shape(0), buf.shape(1)});
    auto res_buf = result.mutable_unchecked<2>();

    for (ssize_t i = 0; i < buf.shape(0); i++)
        for (ssize_t j = 0; j < buf.shape(1); j++)
            res_buf(i, j) = buf(i, j) * buf(i, j);

    return result;
}

PYBIND11_MODULE(rspt_module, m) {
    m.doc() = "2D squaring functions (Python list & NumPy)";

    m.def("square_inplace_list", &square_inplace_list, "Square elements of 2D Python list inplace");
    m.def("square_inplace_numpy", &square_inplace_numpy, "Square elements of 2D NumPy array inplace");
    m.def("square_return_list", &square_return_list, "Return new 2D Python list with squared elements");
    m.def("square_return_numpy", &square_return_numpy, "Return new 2D NumPy array with squared elements");
}
