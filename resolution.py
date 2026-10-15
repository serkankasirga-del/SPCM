"""Resolution of a map from a sharp feature (an edge) in the reflection image.

Draw a line across a sharp edge (e.g. a gold electrode edge). Along the line:
    intensity profile I(d)  ->  first derivative dI/dd  ->  Gaussian fit
The derivative of an edge blurred by the laser spot is the spot profile, so
the FWHM of the fitted Gaussian (2.355 σ) is the resolution of the image.

Works on any 2-D map with ascending axes (also during a scan: pixels not yet
measured are NaN and are skipped).
"""
import os

import numpy as np

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))   # 2.3548


def profile(img, xs, ys, p0, p1, n=None):
    """Bilinear samples of img (rows = ys, columns = xs, both ascending) along
    the straight line p0 -> p1 (same units as xs, ys).
    -> (distance from p0, values); NaN where the line leaves the map or
    touches an unmeasured pixel. Default: 4 samples per map pixel."""
    img = np.asarray(img, float)
    xs, ys = np.asarray(xs, float), np.asarray(ys, float)
    (x0, y0), (x1, y1) = p0, p1
    L = float(np.hypot(x1 - x0, y1 - y0))
    if n is None:
        px = min(abs(xs[1] - xs[0]) if len(xs) > 1 else 1.0,
                 abs(ys[1] - ys[0]) if len(ys) > 1 else 1.0)
        n = int(np.clip(4 * L / px, 20, 4000)) + 1
    t = np.linspace(0.0, 1.0, n)
    X, Y = x0 + t * (x1 - x0), y0 + t * (y1 - y0)
    fx = np.interp(X, xs, np.arange(len(xs)), left=np.nan, right=np.nan)
    fy = np.interp(Y, ys, np.arange(len(ys)), left=np.nan, right=np.nan)
    if len(xs) == 1:
        fx = np.where(np.isclose(X, xs[0]), 0.0, np.nan)
    if len(ys) == 1:
        fy = np.where(np.isclose(Y, ys[0]), 0.0, np.nan)
    out = np.full(n, np.nan)
    ok = np.isfinite(fx) & np.isfinite(fy)
    if ok.any():
        i0 = np.clip(np.floor(fx[ok]).astype(int), 0, max(len(xs) - 2, 0))
        j0 = np.clip(np.floor(fy[ok]).astype(int), 0, max(len(ys) - 2, 0))
        i1 = np.minimum(i0 + 1, len(xs) - 1)
        j1 = np.minimum(j0 + 1, len(ys) - 1)
        a, b = fx[ok] - i0, fy[ok] - j0
        out[ok] = (img[j0, i0] * (1 - a) * (1 - b) + img[j0, i1] * a * (1 - b)
                   + img[j1, i0] * (1 - a) * b + img[j1, i1] * a * b)
    return t * L, out


def _gauss(x, p):
    a, x0, s, c = p
    return a * np.exp(-0.5 * ((x - x0) / s) ** 2) + c


def fit_gaussian(x, y, max_iter=200):
    """Least-squares fit y = a·exp(−(x−x0)²/2σ²) + c (Levenberg–Marquardt).
    -> dict(a, x0, sigma, c, errors, ...) or None."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        return None
    o = np.argsort(x)
    x, y = x[o], y[o]
    m = max(len(x) // 7, 2)                       # baseline from both ends
    c0 = float(np.median(np.concatenate([y[:m], y[-m:]])))
    k = int(np.argmax(np.abs(y - c0)))
    a0 = float(y[k] - c0)
    if a0 == 0:
        return None
    half = np.abs(y - c0) >= abs(a0) / 2          # points above half maximum
    width = float(np.ptp(x[half])) if half.sum() > 1 else float(np.ptp(x)) / 10
    span = float(x[-1] - x[0])
    s_lo, s_hi = max(np.min(np.diff(x)), 1e-9) / 4, span
    p = np.array([a0, x[k], np.clip(width / FWHM_PER_SIGMA, s_lo, s_hi), c0])
    lam = 1e-3

    def resid(p_):
        return y - _gauss(x, p_)

    r = resid(p)
    cost = float(r @ r)
    for _ in range(max_iter):
        e = np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2)
        J = np.column_stack([e, p[0] * e * (x - p[1]) / p[2] ** 2,
                             p[0] * e * (x - p[1]) ** 2 / p[2] ** 3, np.ones_like(x)])
        A, g = J.T @ J, J.T @ r
        try:
            dp = np.linalg.solve(A + lam * np.diag(np.diag(A) + 1e-30), g)
        except np.linalg.LinAlgError:
            break
        pn = p + dp
        pn[1] = np.clip(pn[1], x[0] - 0.5 * span, x[-1] + 0.5 * span)   # stay sensible
        pn[2] = np.clip(abs(pn[2]), s_lo, s_hi)
        rn = resid(pn)
        cn = float(rn @ rn)
        if cn < cost:
            done = cost - cn < 1e-12 * max(cost, 1e-300)
            p, r, cost, lam = pn, rn, cn, lam / 3
            if done:
                break
        else:
            lam *= 4
            if lam > 1e12:
                break
    dof = max(len(x) - 4, 1)
    e = np.exp(-0.5 * ((x - p[1]) / p[2]) ** 2)
    J = np.column_stack([e, p[0] * e * (x - p[1]) / p[2] ** 2,
                         p[0] * e * (x - p[1]) ** 2 / p[2] ** 3, np.ones_like(x)])
    try:
        cov = np.linalg.inv(J.T @ J) * cost / dof
        err = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        err = np.full(4, np.nan)
    return {"a": float(p[0]), "x0": float(p[1]), "sigma": float(abs(p[2])),
            "c": float(p[3]), "a_err": float(err[0]), "x0_err": float(err[1]),
            "sigma_err": float(err[2]), "rms": float(np.sqrt(cost / len(x)))}


def analyse(img, xs, ys, p0, p1, pixel_um=None):
    """Profile along p0 -> p1, its derivative and the Gaussian fit.
    -> dict with dist, intensity, deriv, fit (or None), fwhm_um, fwhm_err_um,
       message. pixel_um: map step, to warn when the result is step-limited."""
    dist, inten = profile(img, xs, ys, p0, p1)
    ok = np.isfinite(inten)
    res = {"dist": dist, "intensity": inten, "deriv": np.full_like(inten, np.nan),
           "fit": None, "fwhm_um": np.nan, "fwhm_err_um": np.nan,
           "p0": tuple(map(float, p0)), "p1": tuple(map(float, p1)),
           "pixel_um": pixel_um}
    if ok.sum() < 8:
        res["message"] = "Not enough measured points along the line."
        return res
    # derivative as a difference across about one map pixel: between pixels
    # the profile is only interpolated, so closer samples give a staircase
    ds = dist[1] - dist[0]
    k = max(1, int(round(0.5 * pixel_um / ds))) if pixel_um else 1
    der = np.full_like(inten, np.nan)
    der[k:-k] = (inten[2 * k:] - inten[:-2 * k]) / (dist[2 * k:] - dist[:-2 * k])
    res["deriv"] = der
    ok = np.isfinite(der)
    if ok.sum() < 8:
        res["message"] = "Not enough measured points along the line."
        return res
    d = dist[ok]
    fit = fit_gaussian(d, der[ok])
    clear = (fit is not None and d[0] <= fit["x0"] <= d[-1]
             and fit["sigma"] < 0.25 * (d[-1] - d[0])          # a peak, not a slope
             and abs(fit["a"]) > 8 * fit["a_err"]               # well determined
             and abs(fit["a"]) > 3 * fit["rms"])                # above the noise
    if not clear:
        res["message"] = ("No clear edge along the line: draw it across one sharp edge "
                          "(e.g. an electrode edge), with flat parts on both sides.")
        return res
    res["fit"] = fit
    res["fwhm_um"] = FWHM_PER_SIGMA * fit["sigma"]
    res["fwhm_err_um"] = FWHM_PER_SIGMA * fit["sigma_err"]
    msg = (f"FWHM = {res['fwhm_um'] * 1000:.0f} ± {res['fwhm_err_um'] * 1000:.0f} nm "
           f"(edge at {fit['x0']:.2f} µm along the line)")
    if pixel_um and res["fwhm_um"] < 2.5 * pixel_um:
        msg += (f". Only {res['fwhm_um'] / pixel_um:.1f} pixels wide: the step "
                f"({pixel_um * 1000:.0f} nm) limits this; scan the edge with a smaller "
                "step for the true resolution")
    res["message"] = msg
    return res


def export(res, path_base, title=""):
    """<base>.csv (distance, intensity, derivative, fit) and <base>.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(os.path.dirname(os.path.abspath(path_base)), exist_ok=True)
    d, inten, der, fit = res["dist"], res["intensity"], res["deriv"], res["fit"]
    gfit = _gauss(d, [fit["a"], fit["x0"], fit["sigma"], fit["c"]]) if fit else \
        np.full_like(d, np.nan)
    with open(path_base + ".csv", "w", encoding="utf-8") as fh:
        fh.write(f"# {title}\n# line from ({res['p0'][0]:.4f}, {res['p0'][1]:.4f}) to "
                 f"({res['p1'][0]:.4f}, {res['p1'][1]:.4f}) um (from the origin)\n")
        fh.write(f"# {res['message']}\n")
        if fit:
            fh.write(f"# fit dI/dd = a exp(-(d-x0)^2 / 2 sigma^2) + c: a={fit['a']:.6e}, "
                     f"x0={fit['x0']:.5f} um, sigma={fit['sigma']:.5f} um, c={fit['c']:.6e}; "
                     f"FWHM = 2.3548 sigma = {res['fwhm_um']:.5f} "
                     f"+- {res['fwhm_err_um']:.5f} um\n")
        fh.write("distance_um,intensity,derivative_per_um,gauss_fit_per_um\n")
        for row in zip(d, inten, der, gfit):
            fh.write(",".join("" if not np.isfinite(v) else f"{v:.6e}" for v in row) + "\n")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(6.0, 5.6), sharex=True)
    a1.plot(d, inten * 1e3, "o-", ms=2.5, lw=1, color="0.25")
    a1.set_ylabel("Reflection (mV)")
    a1.set_title(title or "Resolution from an edge", fontsize=11)
    a2.plot(d, der * 1e3, "o", ms=2.5, color="#1f5fa8", label="dI/dd")
    if fit:
        a2.plot(d, gfit * 1e3, "-", color="#d62728", lw=1.5,
                label=f"Gaussian, FWHM = {res['fwhm_um'] * 1000:.0f} ± "
                      f"{res['fwhm_err_um'] * 1000:.0f} nm")
        a2.legend(fontsize=9)
    a2.axhline(0, color="0.6", lw=0.8)
    a2.set_ylabel("dI/dd (mV/µm)")
    a2.set_xlabel("Distance along the line (µm)")
    fig.tight_layout()
    fig.savefig(path_base + ".png", dpi=200, facecolor="white")
    plt.close(fig)
    return path_base + ".csv", path_base + ".png"
