/**
 * tensorvote_atv.cu
 *
 * pybind11 wrapper for tira::cuda::tensorvote (ATV mode).
 *
 * Python signature:
 *   tensorvote_atv(field, sigma, w=0, stick=True, plate=True, samples=10)
 *     -> np.ndarray shape (H, W, 4), dtype float32
 *
 * The input `field` must be shape (H, W, 4), dtype float32, C-contiguous.
 * Each pixel stores a flat 2x2 symmetric tensor as [a, b, b, c] (row-major).
 *
 * ATV uses a single sigma (sigma1 in the tira::cuda::tensorvote call); sigma2
 * and power are irrelevant for ATV and are set to dummy values here.
 *
 * The output has the same layout: shape (H, W, 4), dtype float32.
 *
 * NOTE: The 2D ATV implementation in tensorvote.h (atv_window / global_tensorvote_atv)
 * does NOT use Boost. Boost is only needed for the 3D analytic plate path.
 * However, because tensorvote.h unconditionally includes both, Boost headers are
 * still required at compile time. See CMakeLists.txt.
 */

// ---- Standard / CUDA headers ----
#include <cuda_runtime.h>

// ---- tira headers ----
#include <tira/functions/tensorvote.h>

// ---- pybind11 headers ----
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

/**
 * Main Python-callable ATV function.
 *
 * Parameters
 * ----------
 * field   : np.ndarray, shape (H, W, 4), dtype float32, C-contiguous
 * sigma   : float — single isotropic decay sigma
 * w       : int   — window side length; 0 = auto (= 2*floor(3*sigma)+1)
 * stick   : bool  — include stick votes
 * plate   : bool  — include plate votes
 * samples : int   — numerical integration samples for ATV plate (>= 1)
 *
 * Returns
 * -------
 * np.ndarray, shape (H, W, 4), dtype float32
 */
py::array_t<float> tensorvote_atv_py(
    py::array_t<float, py::array::c_style | py::array::forcecast> field,
    float sigma,
    int   w       = 0,
    bool  stick   = true,
    bool  plate   = true,
    int   samples = 10)
{
    // ---- Validate input shape ----
    if (field.ndim() != 3 || field.shape(2) != 4)
        throw std::runtime_error("field must have shape (H, W, 4)");

    const int H = static_cast<int>(field.shape(0));
    const int W = static_cast<int>(field.shape(1));

    if (H <= 0 || W <= 0)
        throw std::runtime_error("field must have positive H and W");

    if (sigma <= 0.0f)
        throw std::runtime_error("sigma must be > 0");

    // ---- Auto window size ----
    if (w <= 0) {
        const int pad = static_cast<int>(3.0f * sigma);
        w = 2 * pad + 1;
        if (w < 3) w = 3;
    }
    if (w % 2 == 0) w += 1;

    // ---- Allocate output ----
    py::array_t<float> result({H, W, 4});
    float* out_ptr = result.mutable_data();
    std::memset(out_ptr, 0, sizeof(float) * H * W * 4);

    // ---- Run CUDA ATV voting ----
    // tira::cuda::tensorvote with is_atv=true:
    //   - sigma1 is used as the ATV sigma.
    //   - sigma2 and power are unused in the ATV path but the function signature
    //     requires them; pass safe dummy values (sigma2=0, power=1).
    tira::cuda::tensorvote<float>(
        out_ptr,
        field.data(),
        sigma,
        0.0f,       // sigma2: unused in ATV
        1,          // power:  unused in ATV
        static_cast<unsigned>(w),
        static_cast<unsigned>(H),
        static_cast<unsigned>(W),
        stick,
        plate,
        static_cast<unsigned>(samples),
        true        // is_atv = true
    );

    return result;
}


PYBIND11_MODULE(tensorvote_atv, m) {
    m.doc() = "CUDA-accelerated 2D Analytical Tensor Voting (ATV) via tira::cuda::tensorvote.";

    m.def("tensorvote_atv",
          &tensorvote_atv_py,
          py::arg("field"),
          py::arg("sigma"),
          py::arg("w")       = 0,
          py::arg("stick")   = true,
          py::arg("plate")   = true,
          py::arg("samples") = 10,
          R"doc(
Run CUDA ATV tensor voting on a 2D tensor field.

Parameters
----------
field   : ndarray, shape (H, W, 4), float32
          Flat symmetric tensor field: [T00, T01, T01, T11] per pixel.
sigma   : float — isotropic decay sigma (must be > 0)
w       : int   — window side length (0 = auto: 2*floor(3*sigma)+1)
stick   : bool  — include stick votes
plate   : bool  — include plate votes
samples : int   — numerical integration samples for ATV plate

Returns
-------
ndarray, shape (H, W, 4), float32 — voted tensor field.
)doc");
}
