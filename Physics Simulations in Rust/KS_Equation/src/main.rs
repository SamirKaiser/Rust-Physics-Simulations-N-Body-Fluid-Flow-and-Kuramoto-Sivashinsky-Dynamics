use num_complex::Complex;
use plotters::prelude::*;
use rustfft::{num_traits::Zero, FftPlanner};
use std::f64::consts::PI;

const N: usize = 128;
const L: f64 = 32.0;
const DT: f64 = 0.05;
const NSTEPS: usize = 4000;
const SAVE_EVERY: usize = 4;
const NU: f64 = 0.1;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(N);
    let ifft = planner.plan_fft_inverse(N);

    let dx = L / N as f64;
    let x: Vec<f64> = (0..N).map(|i| i as f64 * dx).collect();

    let mut k = vec![0.0; N];
    for i in 0..N {
        let val = if i <= N / 2 { i as f64 } else { (i as f64) - (N as f64) };
        k[i] = val * (2.0 * PI / L);
    }

    // Initial condition (traveling wave)
    let mut u = vec![0.0; N];
    for i in 0..N {
        let xi = x[i];
        u[i] = 0.3 * (2.0 * PI * xi / L).sin()
            + 0.1 * (4.0 * PI * xi / L).cos()
            + 0.05 * (6.0 * PI * xi / L).sin();
    }

    let mut u_hat: Vec<Complex<f64>> = u.iter().map(|&v| Complex::new(v, 0.0)).collect();
    fft.process(&mut u_hat);

    let mut lk = vec![0.0; N];
    for i in 0..N {
        let k2 = k[i] * k[i];
        lk[i] = -k2 - NU * k2 * k2;
    }

    let e: Vec<Complex<f64>> = lk.iter().map(|&v| Complex::new((DT * v).exp(), 0.0)).collect();
    let e2: Vec<Complex<f64>> = lk.iter().map(|&v| Complex::new((DT * v / 2.0).exp(), 0.0)).collect();

    // ETDRK4 coefficients
    let m = 16;
    let r: Vec<Complex<f64>> = (0..m)
        .map(|mm| Complex::from_polar(1.0, PI * (mm as f64 + 0.5) / m as f64))
        .collect();

    let mut q: Vec<Complex<f64>> = vec![Complex::zero(); N];
    let mut f1: Vec<Complex<f64>> = vec![Complex::zero(); N];
    let mut f2: Vec<Complex<f64>> = vec![Complex::zero(); N];
    let mut f3: Vec<Complex<f64>> = vec![Complex::zero(); N];

    for i in 0..N {
        let mut l_dt_r = vec![Complex::zero(); m];
        for mm in 0..m {
            l_dt_r[mm] = Complex::new(DT * lk[i], 0.0) + r[mm];
        }

        let inv: Vec<Complex<f64>> = l_dt_r.iter().map(|&z| Complex::new(1.0, 0.0) / z).collect();
        let inv2: Vec<Complex<f64>> = l_dt_r.iter().map(|&z| Complex::new(1.0, 0.0) / (z * z)).collect();
        let inv3: Vec<Complex<f64>> = l_dt_r.iter().map(|&z| Complex::new(1.0, 0.0) / (z * z * z)).collect();

        let sum_inv: Complex<f64> = inv.iter().sum();
        let sum_inv2: Complex<f64> = inv2.iter().sum();
        let sum_inv3: Complex<f64> = inv3.iter().sum();

        let dt_c = Complex::new(DT, 0.0);
        let m_c = Complex::new(m as f64, 0.0);

        q[i] = dt_c * (sum_inv / m_c);
        f1[i] = dt_c * ((sum_inv - Complex::new(3.0, 0.0) * sum_inv2 + Complex::new(4.0, 0.0) * sum_inv3) / m_c);
        f2[i] = dt_c * ((Complex::new(2.0, 0.0) * sum_inv - Complex::new(4.0, 0.0) * sum_inv2 + Complex::new(2.0, 0.0) * sum_inv3) / m_c);
        f3[i] = dt_c * ((sum_inv - sum_inv2) / m_c);
    }

    // Time evolution
    let nt_save = NSTEPS / SAVE_EVERY;
    let mut data = vec![vec![0.0; N]; nt_save];

    for step in 0..NSTEPS {
        let n1 = nonlinear(&u_hat, &fft, &ifft, &k);

        let mut a = vec![Complex::zero(); N];
        for i in 0..N {
            a[i] = e2[i] * u_hat[i] + q[i] * n1[i];
        }
        let n2 = nonlinear(&a, &fft, &ifft, &k);

        let mut b = vec![Complex::zero(); N];
        for i in 0..N {
            b[i] = e2[i] * u_hat[i] + q[i] * n2[i];
        }
        let n3 = nonlinear(&b, &fft, &ifft, &k);

        let mut c = vec![Complex::zero(); N];
        for i in 0..N {
            c[i] = e2[i] * a[i] + q[i] * (Complex::new(2.0, 0.0) * n3[i]);
        }
        let n4 = nonlinear(&c, &fft, &ifft, &k);

        for i in 0..N {
            u_hat[i] = e[i] * u_hat[i]
                + f1[i] * n1[i]
                + f2[i] * (n2[i] + n3[i])
                + f3[i] * n4[i];
        }

        if step % SAVE_EVERY == 0 {
            let mut up = u_hat.clone();
            ifft.process(&mut up);
            for i in 0..N {
                data[step / SAVE_EVERY][i] = up[i].re / N as f64;
            }
        }
    }

    println!("✅ Integration complete. Drawing high-res plot...");

    // === Plot settings ===
    let width = 1024;
    let height = 512;
    let root = BitMapBackend::new("ks_spacetime.png", (width, height)).into_drawing_area();
    root.fill(&WHITE)?; // ✅ White background

    let (nt, nx) = (data.len(), N);
    let u_min = data.iter().flatten().fold(f64::INFINITY, |a, &b| a.min(b));
    let u_max = data.iter().flatten().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    let cell_w = width as f64 / nx as f64;
    let cell_h = height as f64 / nt as f64;

    for (ti, row) in data.iter().enumerate() {
        for (xi, &val) in row.iter().enumerate() {
            let norm = (val - u_min) / (u_max - u_min + 1e-12);
            let rgba = HSLColor(norm, 1.0, 0.45).to_rgba(); // ✅ more contrast
            let color = RGBColor(rgba.0, rgba.1, rgba.2);

            let x0 = (xi as f64 * cell_w) as i32;
            let y0 = (ti as f64 * cell_h) as i32;
            let x1 = ((xi + 1) as f64 * cell_w) as i32;
            let y1 = ((ti + 1) as f64 * cell_h) as i32;

            root.draw(&Rectangle::new(
                [(x0, y0), (x1, y1)],
                color.filled(),
            ))?;
        }
    }

    root.present()?;
    println!("✅ Saved high-quality plot: ks_spacetime.png");
    Ok(())
}

// === Nonlinear term ===
fn nonlinear(
    u_hat: &[Complex<f64>],
    fft: &std::sync::Arc<dyn rustfft::Fft<f64>>,
    ifft: &std::sync::Arc<dyn rustfft::Fft<f64>>,
    k: &[f64],
) -> Vec<Complex<f64>> {
    let n = u_hat.len();
    let mut u_phys = u_hat.to_vec();
    ifft.process(&mut u_phys);
    for i in 0..n {
        u_phys[i] /= n as f64;
    }

    let u_sq: Vec<f64> = u_phys.iter().map(|z| z.re * z.re).collect();
    let mut res: Vec<Complex<f64>> = u_sq.iter().map(|&x| Complex::new(x, 0.0)).collect();

    fft.process(&mut res);
    for i in 0..n {
        res[i] *= Complex::new(0.0, k[i]);
        res[i] *= Complex::new(-0.5, 0.0);
    }

    res
}
