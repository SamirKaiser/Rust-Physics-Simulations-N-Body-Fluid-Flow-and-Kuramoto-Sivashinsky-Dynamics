use eframe::egui;
use rayon::prelude::*;

const N: usize = 220;     // grid resolution (wider domain)
const DIFF: f32 = 0.000005;
const VISC: f32 = 0.00002; // lower viscosity → longer vortex street
const DT: f32 = 0.1;

fn ix(x: usize, y: usize) -> usize {
    x + y * N
}

#[derive(Clone)]
struct Fluid {
    size: usize,
    dt: f32,
    diff: f32,
    visc: f32,
    vx: Vec<f32>,
    vy: Vec<f32>,
    vx0: Vec<f32>,
    vy0: Vec<f32>,
    obstacle: Vec<bool>,
    time: f32,
}

impl Fluid {
    fn new(size: usize, diff: f32, visc: f32, dt: f32) -> Self {
        let n2 = size * size;
        let mut obstacle = vec![false; n2];

        // Cylinder near left side
        let center = (size as f32 * 0.25, size as f32 * 0.5);
        let radius = size as f32 * 0.08;
        for j in 0..size {
            for i in 0..size {
                let dx = i as f32 - center.0;
                let dy = j as f32 - center.1;
                if (dx * dx + dy * dy).sqrt() < radius {
                    obstacle[ix(i, j)] = true;
                }
            }
        }

        // Start with a gentle rightward flow
        let vx = vec![1.0; n2];
        let vy = vec![0.0; n2];

        Self {
            size,
            dt,
            diff,
            visc,
            vx: vx.clone(),
            vy,
            vx0: vx.clone(),
            vy0: vec![0.0; n2],
            obstacle,
            time: 0.0,
        }
    }

    fn step(&mut self) {
        let visc = self.visc;
        let diff = self.diff;
        let dt = self.dt;
        let n = self.size;

        diffuse(1, &mut self.vx0, &self.vx, visc, dt, n);
        diffuse(2, &mut self.vy0, &self.vy, visc, dt, n);

        project(&mut self.vx0, &mut self.vy0, &mut self.vx, &mut self.vy, n);

        advect(1, &mut self.vx, &self.vx0, &self.vx0, &self.vy0, dt, n);
        advect(2, &mut self.vy, &self.vy0, &self.vx0, &self.vy0, dt, n);

        project(&mut self.vx, &mut self.vy, &mut self.vx0, &mut self.vy0, n);

        self.time += dt;
        self.apply_inflow_and_obstacle();
    }

    fn apply_inflow_and_obstacle(&mut self) {
        let n = self.size;

        // Smooth inflow ramp
        let inflow_strength = 1.5 / (1.0 + (-0.08 * (self.time - 10.0)).exp());
        for j in 1..n - 1 {
            self.vx[ix(1, j)] = inflow_strength;
            self.vy[ix(1, j)] = 0.0;
        }

        // No-slip obstacle (cylinder)
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                if self.obstacle[ix(i, j)] {
                    self.vx[ix(i, j)] = 0.0;
                    self.vy[ix(i, j)] = 0.0;
                }
            }
        }

        // Outflow (right): open
        for j in 1..n - 1 {
            self.vx[ix(n - 1, j)] = self.vx[ix(n - 2, j)];
            self.vy[ix(n - 1, j)] = self.vy[ix(n - 2, j)];
        }

        // Top/bottom: free-slip
        for i in 0..n {
            self.vx[ix(i, 0)] = self.vx[ix(i, 1)];
            self.vy[ix(i, 0)] = self.vy[ix(i, 1)];
            self.vx[ix(i, n - 1)] = self.vx[ix(i, n - 2)];
            self.vy[ix(i, n - 1)] = self.vy[ix(i, n - 2)];
        }
    }

    fn vorticity(&self) -> Vec<f32> {
        let n = self.size;
        let mut vort = vec![0.0; n * n];
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let dw_dy = (self.vx[ix(i, j + 1)] - self.vx[ix(i, j - 1)]) * 0.5;
                let du_dx = (self.vy[ix(i + 1, j)] - self.vy[ix(i - 1, j)]) * 0.5;
                vort[ix(i, j)] = du_dx - dw_dy;
            }
        }
        vort
    }
}

// ==== Core Solver ====

fn set_bnd(b: usize, x: &mut [f32], n: usize) {
    for i in 1..n - 1 {
        x[ix(i, 0)] = x[ix(i, 1)];
        x[ix(i, n - 1)] = x[ix(i, n - 2)];
        x[ix(0, i)] = x[ix(1, i)];
        x[ix(n - 1, i)] = x[ix(n - 2, i)];
    }
}

fn lin_solve(b: usize, x: &mut [f32], x0: &[f32], a: f32, c: f32, n: usize) {
    for _ in 0..20 {
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                x[ix(i, j)] = (x0[ix(i, j)]
                    + a * (x[ix(i + 1, j)] + x[ix(i - 1, j)] + x[ix(i, j + 1)] + x[ix(i, j - 1)]))
                    / c;
            }
        }
        set_bnd(b, x, n);
    }
}

fn diffuse(b: usize, x: &mut [f32], x0: &[f32], diff: f32, dt: f32, n: usize) {
    let a = dt * diff * (n as f32 - 2.0) * (n as f32 - 2.0);
    lin_solve(b, x, x0, a, 1.0 + 4.0 * a, n);
}

fn advect(b: usize, d: &mut [f32], d0: &[f32], veloc_x: &[f32], veloc_y: &[f32], dt: f32, n: usize) {
    let dt0 = dt * (n as f32 - 2.0);
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let mut x = i as f32 - dt0 * veloc_x[ix(i, j)];
            let mut y = j as f32 - dt0 * veloc_y[ix(i, j)];
            if x < 0.5 { x = 0.5; }
            if x > n as f32 - 1.5 { x = n as f32 - 1.5; }
            if y < 0.5 { y = 0.5; }
            if y > n as f32 - 1.5 { y = n as f32 - 1.5; }
            let i0 = x.floor() as usize;
            let i1 = i0 + 1;
            let j0 = y.floor() as usize;
            let j1 = j0 + 1;
            let s1 = x - i0 as f32;
            let s0 = 1.0 - s1;
            let t1 = y - j0 as f32;
            let t0 = 1.0 - t1;
            d[ix(i, j)] = s0 * (t0 * d0[ix(i0, j0)] + t1 * d0[ix(i0, j1)])
                + s1 * (t0 * d0[ix(i1, j0)] + t1 * d0[ix(i1, j1)]);
        }
    }
    set_bnd(b, d, n);
}

fn project(veloc_x: &mut [f32], veloc_y: &mut [f32], p: &mut [f32], div: &mut [f32], n: usize) {
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            div[ix(i, j)] = -0.5
                * (veloc_x[ix(i + 1, j)] - veloc_x[ix(i - 1, j)] + veloc_y[ix(i, j + 1)]
                    - veloc_y[ix(i, j - 1)])
                / n as f32;
            p[ix(i, j)] = 0.0;
        }
    }
    set_bnd(0, div, n);
    set_bnd(0, p, n);
    lin_solve(0, p, div, 1.0, 4.0, n);
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            veloc_x[ix(i, j)] -= 0.5 * (p[ix(i + 1, j)] - p[ix(i - 1, j)]) * n as f32;
            veloc_y[ix(i, j)] -= 0.5 * (p[ix(i, j + 1)] - p[ix(i, j - 1)]) * n as f32;
        }
    }
    set_bnd(1, veloc_x, n);
    set_bnd(2, veloc_y, n);
}

// ==== Visualization ====

struct VortexApp {
    fluid: Fluid,
}

impl Default for VortexApp {
    fn default() -> Self {
        Self {
            fluid: Fluid::new(N, DIFF, VISC, DT),
        }
    }
}

impl eframe::App for VortexApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.fluid.step();

        egui::CentralPanel::default().show(ctx, |ui| {
            let (response, painter) =
                ui.allocate_painter(ui.available_size(), egui::Sense::hover());
            let rect = response.rect;
            let cell_w = rect.width() / N as f32;
            let cell_h = rect.height() / N as f32;

            let vort = self.fluid.vorticity();
            for j in 0..N {
                for i in 0..N {
                    let v = vort[ix(i, j)];
                    let intensity = (v * 12.0).clamp(-1.0, 1.0); // 🔥 stronger color mapping
                    let color = if self.fluid.obstacle[ix(i, j)] {
                        egui::Color32::DARK_GRAY
                    } else if intensity > 0.0 {
                        egui::Color32::from_rgb((255.0 * intensity) as u8, 40, 40)
                    } else {
                        egui::Color32::from_rgb(40, 40, (-255.0 * intensity) as u8)
                    };
                    let x = rect.left() + i as f32 * cell_w;
                    let y = rect.top() + j as f32 * cell_h;
                    painter.rect_filled(
                        egui::Rect::from_min_size(
                            egui::pos2(x, y),
                            egui::vec2(cell_w, cell_h),
                        ),
                        0.0,
                        color,
                    );
                }
            }
        });

        ctx.request_repaint();
    }
}

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1000.0, 600.0])
            .with_title("Kármán Vortex Street (2D Navier–Stokes)"),
        ..Default::default()
    };

    eframe::run_native(
        "Vortex Street",
        options,
        Box::new(|_cc| Box::<VortexApp>::default()),
    )
}
