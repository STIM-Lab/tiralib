/**
 * tensorvote_tk.cu
 *
 * pybind11 wrapper for tira::cuda::tensorvote (TK mode).
 *
 * Python signature:
 *   tensorvote_tk(field, sigma1, sigma2, power, w=0, stick=True, plate=True, samples=20)
 *     -> np.ndarray shape (H, W, 4), dtype float32
 *
 * The input `field` must be shape (H, W, 4), dtype float32, C-contiguous.
 * Each pixel stores a flat 2x2 symmetric tensor as [a, b, b, c] (row-major).
 * This matches what tira::cuda::tensorvote expects: t_in[idx+0..3].
 *
 * The output has the same layout: shape (H, W, 4), dtype float32.
 *
 * Build: see CMakeLists.txt in this directory.
 *
 * NOTE: This file must be compiled with nvcc (not g++) because it includes
 * tira/functions/tensorvote.h which contains __global__ kernels inside
 * #ifdef __CUDACC__ guards.  The tira:: cuda:: namespace is only active
 * when __CUDACC__ is defined.
 */

// ---- Standard / CUDA headers ----
#include <cuda_runtime.h>

// ---- tira headers ----
// tensorvote.h pulls in tensor.h, eigen.h, callable.h, error.h (all header-only).
// It also pulls in boost headers for the 3D plate analytic form; those are only
// needed for cpu::tk_platevote (3D). The 2D TK path does NOT use Boost at runtime,
// but the include guard does not segregate the 3D section, so Boost headers are
// required at compile time even for the 2D wrapper. See the CMakeLists.txt for
// how to point -I at Boost.
#include <tira/functions/tensorvote.h>

// ---- pybind11 headers ----
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/**
 * Main Python-callable function.
 *
 * Parameters
 * ----------
 * field   : np.ndarray, shape (H, W, 4), dtype float32, C-contiguous
 *           Input tensor field. Each pixel is [T00, T01, T01, T11].
 * sigma1  : float — orthogonal (lateral) sigma
 * sigma2  : float — parallel (axial) sigma
 * power   : int   — refinement exponent (>= 1)
 * w       : int   — voting window side length; 0 = auto (= 2*floor(3*max(s1,s2))+1)
 * stick   : bool  — include stick votes
 * plate   : bool  — include plate votes
 * samples : int   — numerical-integration samples for plate vote (>= 1 for numerical,
 *                   but note: the CUDA kernel ignores plate when samples==0)
 *
 * Returns
 * -------
 * np.ndarray, shape (H, W, 4), dtype float32, C-contiguous
 */
py::array_t<float> tensorvote_tk_py(
    py::array_t<float, py::array::c_style | py::array::forcecast> field,
    float sigma1,
    float sigma2,
    int   power,
    int   w       = 0,
    bool  stick   = true,
    bool  plate   = true,
    int   samples = 20)
{
    // ---- Validate input shape ----
    if (field.ndim() != 3 || field.shape(2) != 4)
        throw std::runtime_error("field must have shape (H, W, 4)");

    const int H = static_cast<int>(field.shape(0));
    const int W = static_cast<int>(field.shape(1));

    if (H <= 0 || W <= 0)
        throw std::runtime_error("field must have positive H and W");

    if (sigma1 <= 0.0f)
        throw std::runtime_error("sigma1 must be > 0");

    if (power < 1) power = 1;

    // ---- Auto window size ----
    if (w <= 0) {
        const float sigmax = (sigma2 > sigma1) ? sigma2 : sigma1;
        const int pad = static_cast<int>(3.0f * sigmax);
        w = 2 * pad + 1;
        if (w < 3) w = 3;
    }
    // Window must be odd
    if (w % 2 == 0) w += 1;

    // ---- Allocate output (zero-initialized, host) ----
    py::array_t<float> result({H, W, 4});
    float* out_ptr = result.mutable_data();
    std::memset(out_ptr, 0, sizeof(float) * H * W * 4);

    // ---- Run CUDA TK voting ----
    // tira::cuda::tensorvote signature (2D, is_atv=false):
    //   tensorvote<float>(t_out, t_in, sigma1, sigma2, power, w, shape0, shape1,
    //                     stick, plate, samples, is_atv=false)
    // Both t_out and t_in are host pointers here; the function copies to/from device
    // internally via _dev_readonly / _dev_writeable helpers.
    tira::cuda::tensorvote<float>(
        out_ptr,
        field.data(),
        sigma1,
        sigma2,
        power,
        static_cast<unsigned>(w),
        static_cast<unsigned>(H),
        static_cast<unsigned>(W),
        stick,
        plate,
        static_cast<unsigned>(samples),
        false   // is_atv = false -> TK mode
    );

    return result;
}


/**
 * tensorvote_tk_prebuilt_py
 *
 * Like tensorvote_tk_py but accepts a pre-built neighbor table from Python
 * (constructed via numpy) instead of calling cpu::build_neighbors2d internally.
 * This eliminates the O(w^2 * N) scalar C++ table-build loop that dominates
 * at large sigma or power.
 *
 * Parameters
 * ----------
 * field    : ndarray, shape (H, W, 4), float32, C-contiguous
 * nb_buf   : buffer whose raw bytes hold nb_count tira::Neighbor2D structs
 *            (48 bytes each: int32 du, int32 dv, float32 dx/dy/c1/c2/Ra/Rb/Rc/Va/Vb/Vc)
 * nb_count : number of entries in nb_buf
 * w        : voting window side length (must be odd, > 0)
 * sigma1   : float — orthogonal sigma (must be > 0)
 * sigma2   : float — axial sigma
 * power    : int   — refinement exponent
 * stick    : bool
 * plate    : bool
 * samples  : int   — plate numerical integration samples
 */
py::array_t<float> tensorvote_tk_prebuilt_py(
    py::array_t<float, py::array::c_style | py::array::forcecast> field,
    py::buffer   nb_buf,
    int          nb_count,
    int          w,
    float        sigma1,
    float        sigma2  = 0.0f,
    int          power   = 1,
    bool         stick   = true,
    bool         plate   = true,
    int          samples = 20)
{
    if (field.ndim() != 3 || field.shape(2) != 4)
        throw std::runtime_error("field must have shape (H, W, 4)");

    const int H = static_cast<int>(field.shape(0));
    const int W = static_cast<int>(field.shape(1));
    if (H <= 0 || W <= 0) throw std::runtime_error("field must have positive H and W");
    if (sigma1 <= 0.0f) throw std::runtime_error("sigma1 must be > 0");
    if (power < 1) power = 1;
    if (w <= 0 || w % 2 == 0) throw std::runtime_error("w must be a positive odd integer");

    py::buffer_info nb_info = nb_buf.request();
    const std::size_t nb_bytes_expected =
        static_cast<std::size_t>(nb_count) * sizeof(tira::Neighbor2D);
    const std::size_t nb_bytes_actual =
        static_cast<std::size_t>(nb_info.itemsize) * static_cast<std::size_t>(nb_info.size);
    if (nb_bytes_actual != nb_bytes_expected)
        throw std::runtime_error(
            "nb_buf byte size does not match nb_count * sizeof(Neighbor2D) (48)");

    const tira::Neighbor2D* host_nb =
        reinterpret_cast<const tira::Neighbor2D*>(nb_info.ptr);

    py::array_t<float> result({H, W, 4});
    float* out_ptr = result.mutable_data();
    std::memset(out_ptr, 0, sizeof(float) * H * W * 4);

    tira::cuda::tensorvote_prebuilt<float>(
        out_ptr,
        field.data(),
        host_nb,
        nb_count,
        sigma1,
        sigma2,
        power,
        static_cast<unsigned>(w),
        static_cast<unsigned>(H),
        static_cast<unsigned>(W),
        stick,
        plate,
        static_cast<unsigned>(samples));

    return result;
}


PYBIND11_MODULE(tensorvote_tk, m) {
    m.doc() = "CUDA-accelerated 2D Tensor-Kernel (TK) tensor voting via tira::cuda::tensorvote.";

    m.def("tensorvote_tk",
          &tensorvote_tk_py,
          py::arg("field"),
          py::arg("sigma1"),
          py::arg("sigma2")  = 0.0f,
          py::arg("power")   = 1,
          py::arg("w")       = 0,
          py::arg("stick")   = true,
          py::arg("plate")   = true,
          py::arg("samples") = 20,
          R"doc(
Run CUDA TK tensor voting on a 2D tensor field.

Parameters
----------
field   : ndarray, shape (H, W, 4), float32
          Flat symmetric tensor field. Each pixel stores [T00, T01, T01, T11].
sigma1  : float — orthogonal (lateral) decay sigma (must be > 0)
sigma2  : float — parallel (axial) decay sigma (0 = no parallel channel)
power   : int   — refinement exponent >= 1
w       : int   — window side length (0 = auto: 2*floor(3*max(sigma1,sigma2))+1)
stick   : bool  — include stick votes
plate   : bool  — include plate votes
samples : int   — numerical integration samples for plate (>= 1)

Returns
-------
ndarray, shape (H, W, 4), float32 — voted tensor field, same layout as input.
)doc");

    m.def("tensorvote_tk_prebuilt",
          &tensorvote_tk_prebuilt_py,
          py::arg("field"),
          py::arg("nb_buf"),
          py::arg("nb_count"),
          py::arg("w"),
          py::arg("sigma1"),
          py::arg("sigma2")  = 0.0f,
          py::arg("power")   = 1,
          py::arg("stick")   = true,
          py::arg("plate")   = true,
          py::arg("samples") = 20,
          "TK voting with a pre-built neighbor table (48 bytes/entry). "
          "Skips the internal cpu::build_neighbors2d call; use _build_nb_numpy() "
          "from cuda_wrappers.py to build the table with vectorized numpy.");
}
